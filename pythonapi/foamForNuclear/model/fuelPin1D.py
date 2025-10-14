import foamForNuclear.boundaryConditions as bc
from foamForNuclear.field import Dimension, Field, GapGas
from foamForNuclear.mesh.blockMesh import BlockMesh, Face
from foamForNuclear.offbeat import OffbeatSolver
from foamForNuclear.timeFolder import TimeFolder
from .model import Model


class FuelPin1D(Model):
    """
    Model of a 1D fuel pin using the OFFBEAT solver.
    """
    def __init__(
            self,
            rInnerFuel: float,
            rOuterFuel: float,
            rInnerClad: float,
            rOuterClad: float,
            lFuel: float,
            lPlenum: float,
            wedgeAngle: float,
            nrFuel: int,
            nrClad: int,
            nzFuel: int,
            nzPlenum: int,
            zFuelBottom: float=0,
            scale: float=1
        ):
        super().__init__()

        self.rInnerFuel = rInnerFuel
        self.rOuterFuel = rOuterFuel
        self.rInnerClad = rInnerClad
        self.rOuterClad = rOuterClad
        self.lFuel = lFuel
        self.lPlenum = lPlenum
        self.wedgeAngle = wedgeAngle
        self.nrFuel = nrFuel
        self.nrClad = nrClad
        self.nzFuel = nzFuel
        self.nzPlenum = nzPlenum
        self.zFuelBottom = zFuelBottom

        self.scale = scale

        self.settings.application = 'offbeat'
        self.mesh = None

        self.rescale(scale)

        self.create_mesh()
        self.create_time_folder()

        self.create_solver()


    def create_mesh(self):
        self.mesh = BlockMesh()

        # Create internal meshes
        fuel = self.mesh.createWedge(
            'fuel',
            self.rInnerFuel, self.rOuterFuel,
            self.zFuelBottom, self.zFuelBottom + self.lFuel,
            self.wedgeAngle,
            nr=self.nrFuel,
            nz=self.nzFuel
        )

        cladding = self.mesh.createWedge(
            'cladding',
            self.rInnerClad, self.rOuterClad,
            self.zFuelBottom, self.zFuelBottom + self.lFuel,
            self.wedgeAngle,
            nr=self.nrClad,
            nz=self.nzFuel
        )

        claddingPlenum = self.mesh.extrudeTop(
            [cladding], 'cladding', self.lPlenum, self.nzPlenum
        )


        # Create boundary conditions for fuel
        fuelBottom = Face("fuelBottom", boundaryType="empty")
        fuelBottom.addSubFace(fuel.bottomFace())

        fuelTop = Face("fuelTop", boundaryType="empty")
        fuelTop.addSubFace(fuel.topFace())

        fuelFront = Face("fuelFront", boundaryType="wedge")
        fuelFront.addSubFace(fuel.frontFace())

        fuelBack = Face("fuelBack", boundaryType="wedge")
        fuelBack.addSubFace(fuel.backFace())

        fuelOuter = Face("fuelOuter", boundaryType="regionCoupledOFFBEAT")
        fuelOuter.addSubFace(fuel.rightFace())
        fuelOuter.extraParameters = {
            "neighbourPatch": "cladInner",
            "neighbourRegion": "region0",
            "owner": "true",
            "updateAMI": "true"
        }

        self.mesh.addBoundary(fuelBottom)
        self.mesh.addBoundary(fuelTop)
        self.mesh.addBoundary(fuelFront)
        self.mesh.addBoundary(fuelBack)
        self.mesh.addBoundary(fuelOuter)

        if (self.rInnerFuel > 0):
            fuelInner = Face("fuelInner")
            fuelInner.addSubFace(fuel.leftFace())
            self.mesh.addBoundary(fuelInner)


        # Create boundary conditions for cladding
        cladBottom = Face("cladBottom", boundaryType="empty")
        cladBottom.addSubFace(cladding.bottomFace())

        cladTop = Face("cladTop", boundaryType="empty")
        cladTop.addSubFace(claddingPlenum.topFace())

        cladFront = Face("cladFront", boundaryType="wedge")
        cladFront.addSubFace(cladding.frontFace())
        cladFront.addSubFace(claddingPlenum.frontFace())

        cladBack = Face("cladBack", boundaryType="wedge")
        cladBack.addSubFace(cladding.backFace())
        cladBack.addSubFace(claddingPlenum.backFace())

        cladInner = Face("cladInner", boundaryType="regionCoupledOFFBEAT")
        cladInner.addSubFace(cladding.leftFace())
        cladInner.addSubFace(claddingPlenum.leftFace())
        cladInner.extraParameters = {
            "neighbourPatch": "fuelOuter",
            "neighbourRegion": "region0",
            "owner": "false",
            "updateAMI": "true"
        }

        cladOuter = Face("cladOuter")
        cladOuter.addSubFace(cladding.rightFace())
        cladOuter.addSubFace(claddingPlenum.rightFace())

        self.mesh.addBoundary(cladBottom)
        self.mesh.addBoundary(cladTop)
        self.mesh.addBoundary(cladFront)
        self.mesh.addBoundary(cladBack)
        self.mesh.addBoundary(cladInner)
        self.mesh.addBoundary(cladOuter)


    def create_time_folder(self):
        timeFolder0 = TimeFolder(0)
        self.timeFolders.append(timeFolder0)

        # Temperature
        T = Field("T")
        T.dimensions = Dimension(default='T')
        T.internalField = 300
        T.set_boundary_condition('".*Front|.*Back"', bc.Wedge())
        T.set_boundary_condition('fuelOuter', bc.FuelRodGap(
            value=T.internalField,
            roughness=2.2e-6
        ))
        T.set_boundary_condition('cladInner', bc.FuelRodGap(
            value=T.internalField,
            roughness=0.5e-6
        ))
        T.set_boundary_condition('cladOuter', bc.UniformFixedValue(
            uniformValue=[
                (0,        300),
                (3600,     600),
                (3.15E+07, 600),
            ],
        ))
        timeFolder0.append(T)

        # Neutron flux
        neutronFlux0 = Field("neutronFlux0")
        neutronFlux0.dimensions = Dimension(default='neutronFlux')
        neutronFlux0.internalField = 0
        neutronFlux0.set_boundary_condition('".*Front|.*Back"', bc.Wedge())
        neutronFlux0.set_boundary_condition('"cladInner|cladOuter"', bc.FixedValue(0))
        neutronFlux0.set_boundary_condition('fuelOuter', bc.FixedValue(1))
        timeFolder0.append(neutronFlux0)

        # Gap gas
        gapGas = GapGas('fromModel')
        gapGas.gasPressure = 2.25e6
        timeFolder0.append(gapGas)


    def create_solver(self):
        mechSolver = OffbeatSolver(
            solver="offbeat",
            mesh=self.mesh
        )

        self.solvers.append(mechSolver)


    def rescale(self, scale: float):
        self.scale = scale
        self.rInnerFuel *= scale
        self.rOuterFuel *= scale
        self.rInnerClad *= scale
        self.rOuterClad *= scale
        self.lFuel *= scale
        self.lPlenum *= scale
