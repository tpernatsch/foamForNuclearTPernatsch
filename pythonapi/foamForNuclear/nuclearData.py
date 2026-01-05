import foamlib
import numpy as np
from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.timeProfile import TimeProfile


_XS_VARIABLE_LAW_TYPES = {'lin', 'log', 'sqrt'}
_XS_TYPES = {
    'sigmaRemoval', 'nuSigmaEff', 'sigmaPow', 'scatteringMatrixP0',
    'scatteringMatrixP1', 'scatteringMatrixP2', 'scatteringMatrixP3',
    'scatteringMatrixP4', 'scatteringMatrixP5', 'discFactor', 'chiPrompt',
    'chiDelayed', 'IV', 'lambda', 'Beta', 'D', 'integralFlux'
}
_XS_UNITS_DICT = {
    "sigmaRemoval": "1/m",
    "nuSigmaEff": "neutrons/m",
    "sigmaPow": "J/m",
    "scatteringMatrixP0": "1/m",
    "scatteringMatrixP1": "1/m",
    "scatteringMatrixP2": "1/m",
    "scatteringMatrixP3": "1/m",
    "scatteringMatrixP4": "1/m",
    "scatteringMatrixP5": "1/m",
    "discFactor": "-",
    "chiPrompt": "-",
    "chiDelayed": "-",
    "IV": "s/m",
    "lambda": "1/s",
    "Beta": "-",
    "D": "m",
    "integralFlux": "-",
}


class NuclearDataZone(OpenFOAMDict):
    """
    Nuclear data object for a single zone in a specific state.

    Example of use with `read_from_openmc`::

        openmcZone.read_from_openmc(
            outputFilename="../OpenMC/statepoint.500.h5",
            domains=openmc.Cell(name='Fuel'),
            energyGroups=openmc.mgxs.EnergyGroups([0., 0.625, 20.0e6]),
            delayedGroups=list(range(1,7)),
            rootDomain=openmc.Universe(0, name='root universe')
        )

    Example of use with `read_from_serpent`::

        serpentZone.read_from_serpent(
            outputFilename="../Serpent/main_res.m",
            universe='uFuel'
        )

    Parameters
    ----------
    fuelFraction : float
        Volume fraction of fuel in a cell (default `1`). It is possible to use
        Honeycomb https://foam-for-nuclear.gitlab.io/honeycomb/ to compute this
        quantity.
    inverseVelocity : NonUniformList
        Corresponds to the parameter `IV`
    diffusionCoefficient : NonUniformList
        Corresponds to the parameter `D`
    nuFissionXS : NonUniformList
        Corresponds to the parameter `nuSigmaEff`
    powerXS : NonUniformList
        Corresponds to the parameter `sigmaPow` and is equal to
        kappa * Sigma_fission
    scatteringMatrixP0 : Matrix
        Corresponds to the parameter `scatteringMatrixP0`
    scatteringMatrixP1 : Matrix
        Corresponds to the parameter `scatteringMatrixP1`
    scatteringMatrixP2 : Matrix
        Corresponds to the parameter `scatteringMatrixP2`
    scatteringMatrixP3 : Matrix
        Corresponds to the parameter `scatteringMatrixP3`
    scatteringMatrixP4 : Matrix
        Corresponds to the parameter `scatteringMatrixP4`
    scatteringMatrixP5 : Matrix
        Corresponds to the parameter `scatteringMatrixP5`
    removalXS : NonUniformList
        Corresponds to the parameter `sigmaRemoval`
    chiPrompt : NonUniformList
        Corresponds to the parameter `chiPrompt`
    chiDelayed : NonUniformList
        Corresponds to the parameter `chiDelayed`
    delayedFraction : NonUniformList
        Corresponds to the parameter `Beta`
    decayConstant : NonUniformList
        Corresponds to the parameter `lambda`
    discFactor : NonUniformList
        Corresponds to the parameter `discFactor`
    integralFlux : NonUniformList
        Corresponds to the parameter `integralFlux`
    description : str
        Add a description in the `nuclearData` file during export
        (default `None`).
    groupStructure : list[float]
        Optional: Only used for information and plotting. Energy group structure
        in increasing order in MeV (default `None`).
        If `read_from_serpent` or `read_from_openmc` is used, the
        `groupStructure` is automatically filled with the data in the results
        file.

        Example for 2 energy groups::

            groupStructure = [0, 0.625e-6, 20]

    Attributes
    ----------
    energyGroups : int
        Number of energy groups based on the length of `inverseVelocity`.
    precGroups : int
        Number of precursor/delayed neutron groups based on the length of
        `decayConstant`.
    """
    def __init__(
            self,
            name,
            fuelFraction: float=1,
            inverseVelocity=NonUniformList([]),
            diffusionCoefficient=NonUniformList([]),
            nuFissionXS=NonUniformList([]),
            powerXS=NonUniformList([]),
            scatteringMatrixP0=Matrix([]),
            scatteringMatrixP1=Matrix([]),
            scatteringMatrixP2=Matrix([]),
            scatteringMatrixP3=Matrix([]),
            scatteringMatrixP4=Matrix([]),
            scatteringMatrixP5=Matrix([]),
            removalXS=NonUniformList([]),
            chiPrompt=NonUniformList([]),
            chiDelayed=NonUniformList([]),
            delayedFraction=NonUniformList([]),
            decayConstant=NonUniformList([]),
            discFactor=NonUniformList([]),
            integralFlux=NonUniformList([]),
            description: str=None,
            groupStructure: list[float]=None
        ):
        super().__init__(name=name)
        self.description = description
        self.groupStructure = groupStructure
        self.fuelFraction = fuelFraction
        self.inverseVelocity = inverseVelocity
        self.diffusionCoefficient = diffusionCoefficient
        self.nuFissionXS = nuFissionXS
        self.powerXS = powerXS
        self.scatteringMatrixP0 = scatteringMatrixP0
        self.scatteringMatrixP1 = scatteringMatrixP1
        self.scatteringMatrixP2 = scatteringMatrixP2
        self.scatteringMatrixP3 = scatteringMatrixP3
        self.scatteringMatrixP4 = scatteringMatrixP4
        self.scatteringMatrixP5 = scatteringMatrixP5
        self.removalXS = removalXS
        self.chiPrompt = chiPrompt
        self.chiDelayed = chiDelayed
        self.delayedFraction = delayedFraction
        self.decayConstant = decayConstant
        self.discFactor = discFactor
        self.integralFlux = integralFlux


    def __repr__(self, depth: int=0):
        # Check main data
        for data, dataName in [
            (self.inverseVelocity, "inverseVelocity"),
            (self.diffusionCoefficient, "diffusionCoefficient"),
            (self.nuFissionXS, "nuFissionXS"),
            (self.powerXS, "powerXS"),
            (self.removalXS, "removalXS"),
            (self.chiPrompt, "chiPrompt"),
            (self.chiDelayed, "chiDelayed"),
            (self.discFactor, "discFactor"),
            (self.integralFlux, "integralFlux"),
        ]:
            if (len(data) != self.energyGroups):
                msg = f"Length of {dataName} = {len(data)} which is different from energyGroups = {self.energyGroups}"
                raise IndexError(msg)

        # Check delayed neutrons
        for data, dataName in [
            (self.delayedFraction, "delayedFraction"),
            (self.decayConstant, "decayConstant"),
        ]:
            if (len(data) != self.precGroups):
                msg = f"Length of {dataName} = {len(data)} which is different from precGroups = {self.precGroups}"
                raise IndexError(msg)

        # Check matrices
        for data, dataName in [
            (self.scatteringMatrixP0, "scatteringMatrixP0"),
            (self.scatteringMatrixP1, "scatteringMatrixP1"),
            (self.scatteringMatrixP2, "scatteringMatrixP2"),
            (self.scatteringMatrixP3, "scatteringMatrixP3"),
            (self.scatteringMatrixP4, "scatteringMatrixP4"),
            (self.scatteringMatrixP5, "scatteringMatrixP5"),
        ]:
            if (
                not data.is_empty and (
                    len(data) != self.energyGroups or
                    len(data[0]) != self.energyGroups
                )
            ):
                msg = f"Length of {dataName} = ({len(data)}, {len(data[0])}) which is different from energyGroups = ({self.energyGroups}, {self.energyGroups})"
                raise IndexError(msg)

        self.__setitem__("fuelFraction", self.fuelFraction)
        self.__setitem__("IV", self.inverseVelocity.__repr__(noBreak=True))
        self.__setitem__("D", self.diffusionCoefficient.__repr__(noBreak=True))
        self.__setitem__("nuSigmaEff", self.nuFissionXS.__repr__(noBreak=True))
        self.__setitem__("sigmaPow", self.powerXS.__repr__(noBreak=True))
        self.__setitem__("scatteringMatrixP0", self.scatteringMatrixP0.__repr__(depth=depth+1))
        if (not self.scatteringMatrixP1.is_empty):
            self.__setitem__("scatteringMatrixP1", self.scatteringMatrixP1.__repr__(depth=depth+1))
        if (not self.scatteringMatrixP2.is_empty):
            self.__setitem__("scatteringMatrixP2", self.scatteringMatrixP2.__repr__(depth=depth+1))
        if (not self.scatteringMatrixP3.is_empty):
            self.__setitem__("scatteringMatrixP3", self.scatteringMatrixP3.__repr__(depth=depth+1))
        if (not self.scatteringMatrixP4.is_empty):
            self.__setitem__("scatteringMatrixP4", self.scatteringMatrixP4.__repr__(depth=depth+1))
        if (not self.scatteringMatrixP5.is_empty):
            self.__setitem__("scatteringMatrixP5", self.scatteringMatrixP5.__repr__(depth=depth+1))
        self.__setitem__("sigmaRemoval", self.removalXS.__repr__(noBreak=True))
        self.__setitem__("chiPrompt", self.chiPrompt.__repr__(noBreak=True))
        self.__setitem__("chiDelayed", self.chiDelayed.__repr__(noBreak=True))
        self.__setitem__("Beta", self.delayedFraction.__repr__(noBreak=True))
        self.__setitem__("lambda", self.decayConstant.__repr__(noBreak=True))
        self.__setitem__("discFactor", self.discFactor.__repr__(noBreak=True))
        self.__setitem__("integralFlux", self.integralFlux.__repr__(noBreak=True))

        text = ""
        if (self.description is not None):
            text += f"// {self.description}\n{depth*tab}"
        return(text + super().__repr__(depth))

    @property
    def energyGroups(self) -> int:
        return(len(self.inverseVelocity))

    @property
    def precGroups(self) -> int:
        return(len(self.decayConstant))

    @property
    def fuelFraction(self):
        return self._fuelFraction

    @fuelFraction.setter
    def fuelFraction(self, fuelFraction):
        check_type("fuelFraction", fuelFraction, (int, float))
        check_positive("fuelFraction", fuelFraction)
        self._fuelFraction = fuelFraction

    @property
    def inverseVelocity(self):
        return self._inverseVelocity

    @inverseVelocity.setter
    def inverseVelocity(self, inverseVelocity):
        check_type("inverseVelocity", inverseVelocity, (list, List, NonUniformList, np.ndarray))
        if (isinstance(inverseVelocity, (list, List, np.ndarray))):
            self._inverseVelocity = NonUniformList(inverseVelocity)
        else:
            self._inverseVelocity = inverseVelocity

    @property
    def diffusionCoefficient(self):
        return self._diffusionCoefficient

    @diffusionCoefficient.setter
    def diffusionCoefficient(self, diffusionCoefficient):
        check_type("diffusionCoefficient", diffusionCoefficient, (list, List, np.ndarray, NonUniformList))
        if (isinstance(diffusionCoefficient, (list, List, np.ndarray))):
            self._diffusionCoefficient = NonUniformList(diffusionCoefficient)
        else:
            self._diffusionCoefficient = diffusionCoefficient

    @property
    def nuFissionXS(self):
        return self._nuSigmaEff

    @nuFissionXS.setter
    def nuFissionXS(self, nuFissionXS):
        check_type("nuFissionXS", nuFissionXS, (list, List, np.ndarray, NonUniformList))
        if (isinstance(nuFissionXS, (list, List, np.ndarray))):
            self._nuSigmaEff = NonUniformList(nuFissionXS)
        else:
            self._nuSigmaEff = nuFissionXS

    @property
    def powerXS(self):
        return self._sigmaPow

    @powerXS.setter
    def powerXS(self, powerXS):
        check_type("powerXS", powerXS, (list, List, np.ndarray, NonUniformList))
        if (isinstance(powerXS, (list, List, np.ndarray))):
            self._sigmaPow = NonUniformList(powerXS)
        else:
            self._sigmaPow = powerXS

    @property
    def removalXS(self):
        return self._sigmaRemoval

    @removalXS.setter
    def removalXS(self, removalXS):
        check_type("removalXS", removalXS, (list, List, np.ndarray, NonUniformList))
        if (isinstance(removalXS, (list, List, np.ndarray))):
            self._sigmaRemoval = NonUniformList(removalXS)
        else:
            self._sigmaRemoval = removalXS

    @property
    def chiPrompt(self):
        return self._chiPrompt

    @chiPrompt.setter
    def chiPrompt(self, chiPrompt):
        check_type("chiPrompt", chiPrompt, (list, List, np.ndarray, NonUniformList))
        if (isinstance(chiPrompt, (list, List, np.ndarray))):
            self._chiPrompt = NonUniformList(chiPrompt)
        else:
            self._chiPrompt = chiPrompt

    @property
    def chiDelayed(self):
        return self._chiDelayed

    @chiDelayed.setter
    def chiDelayed(self, chiDelayed):
        check_type("chiDelayed", chiDelayed, (list, List, np.ndarray, NonUniformList))
        if (isinstance(chiDelayed, (list, List, np.ndarray))):
            self._chiDelayed = NonUniformList(chiDelayed)
        else:
            self._chiDelayed = chiDelayed

    @property
    def delayedFraction(self):
        return self._Beta

    @delayedFraction.setter
    def delayedFraction(self, delayedFraction):
        check_type("delayedFraction", delayedFraction, (list, List, np.ndarray, NonUniformList))
        if (isinstance(delayedFraction, (list, List, np.ndarray))):
            self._Beta = NonUniformList(delayedFraction)
        else:
            self._Beta = delayedFraction

    @property
    def decayConstant(self):
        return self._decayConstant

    @decayConstant.setter
    def decayConstant(self, decayConstant):
        check_type("decayConstant", decayConstant, (list, List, np.ndarray, NonUniformList))
        if (isinstance(decayConstant, (list, List, np.ndarray))):
            self._decayConstant = NonUniformList(decayConstant)
        else:
            self._decayConstant = decayConstant

    @property
    def discFactor(self):
        return self._discFactor

    @discFactor.setter
    def discFactor(self, discFactor):
        check_type("discFactor", discFactor, (list, List, np.ndarray, NonUniformList))
        if (isinstance(discFactor, (list, List, np.ndarray))):
            self._discFactor = NonUniformList(discFactor)
        else:
            self._discFactor = discFactor

    @property
    def integralFlux(self):
        return self._integralFlux

    @integralFlux.setter
    def integralFlux(self, integralFlux):
        check_type("integralFlux", integralFlux, (list, List, np.ndarray, NonUniformList))
        if (isinstance(integralFlux, (list, List, np.ndarray))):
            self._integralFlux = NonUniformList(integralFlux)
        else:
            self._integralFlux = integralFlux

    @property
    def scatteringMatrixP0(self):
        return self._scatteringMatrixP0

    @scatteringMatrixP0.setter
    def scatteringMatrixP0(self, scatteringMatrixP0):
        check_type("scatteringMatrixP0", scatteringMatrixP0, (list, np.ndarray, Matrix))
        if (isinstance(scatteringMatrixP0, (list, np.ndarray))):
            self._scatteringMatrixP0 = Matrix(scatteringMatrixP0)
        else:
            self._scatteringMatrixP0 = scatteringMatrixP0

    @property
    def scatteringMatrixP1(self):
        return self._scatteringMatrixP1

    @scatteringMatrixP1.setter
    def scatteringMatrixP1(self, scatteringMatrixP1):
        check_type("scatteringMatrixP1", scatteringMatrixP1, (list, np.ndarray, Matrix))
        if (isinstance(scatteringMatrixP1, (list, np.ndarray))):
            self._scatteringMatrixP1 = Matrix(scatteringMatrixP1)
        else:
            self._scatteringMatrixP1 = scatteringMatrixP1

    @property
    def scatteringMatrixP2(self):
        return self._scatteringMatrixP2

    @scatteringMatrixP2.setter
    def scatteringMatrixP2(self, scatteringMatrixP2):
        check_type("scatteringMatrixP2", scatteringMatrixP2, (list, np.ndarray, Matrix))
        if (isinstance(scatteringMatrixP2, (list, np.ndarray))):
            self._scatteringMatrixP2 = Matrix(scatteringMatrixP2)
        else:
            self._scatteringMatrixP2 = scatteringMatrixP2

    @property
    def scatteringMatrixP3(self):
        return self._scatteringMatrixP3

    @scatteringMatrixP3.setter
    def scatteringMatrixP3(self, scatteringMatrixP3):
        check_type("scatteringMatrixP3", scatteringMatrixP3, (list, np.ndarray, Matrix))
        if (isinstance(scatteringMatrixP3, (list, np.ndarray))):
            self._scatteringMatrixP3 = Matrix(scatteringMatrixP3)
        else:
            self._scatteringMatrixP3 = scatteringMatrixP3

    @property
    def scatteringMatrixP4(self):
        return self._scatteringMatrixP4

    @scatteringMatrixP4.setter
    def scatteringMatrixP4(self, scatteringMatrixP4):
        check_type("scatteringMatrixP4", scatteringMatrixP4, (list, np.ndarray, Matrix))
        if (isinstance(scatteringMatrixP4, (list, np.ndarray))):
            self._scatteringMatrixP4 = Matrix(scatteringMatrixP4)
        else:
            self._scatteringMatrixP4 = scatteringMatrixP4

    @property
    def scatteringMatrixP5(self):
        return self._scatteringMatrixP5

    @scatteringMatrixP5.setter
    def scatteringMatrixP5(self, scatteringMatrixP5):
        check_type("scatteringMatrixP5", scatteringMatrixP5, (list, np.ndarray, Matrix))
        if (isinstance(scatteringMatrixP5, (list, np.ndarray))):
            self._scatteringMatrixP5 = Matrix(scatteringMatrixP5)
        else:
            self._scatteringMatrixP5 = scatteringMatrixP5

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, description):
        check_type("description", description, str, none_ok=True)
        self._description = description


    def read_from_serpent(
            self,
            outputFilename: str,
            universe: str='0',
            step: int=0,
            burnup: float=None,
            days: float=None,
            isEffective: bool=True,
            serpentToolsResfile=None
        ) -> None:
        """
        Extract the nuclear data from Serpent `.m` output file for a single
        universe and burnup step.

        Parameters
        ----------
        outputFilename : str
            Serpent output file name
        universe : str
            Name of the Serpent universe to extract the MGXS
        step : int
            Burnup step (default to 0)
        burnup : float
            Burnup [MWd/kgU] of the desired universe
        days : float
            Time [days] of the desired universe
        isEffective : bool
            Flag to extract the effective or analogue values (default to True)
        serpentToolsResfile :
            serpentTools results file already open
        """
        import serpentTools

        # Unit conversion from Serpent to OpenFOAM
        mev2j = 1.602176487e-13
        cm2m = 0.01

        res = serpentToolsResfile
        if (serpentToolsResfile is None):
            res = serpentTools.read(outputFilename)

        # Extract universe
        univ = res.getUniv(universe, burnup=burnup, index=step, timeDays=days)

        burnup = univ.bu
        step = univ.step
        days = univ.day

        # Reverse to be in increasing order
        self.groupStructure = univ.groups[::-1]

        ng = univ.numGroups

        # Extract and convert nuclear data
        inverseVelocity = univ.infExp['infInvv'] / cm2m
        diffcoeff = univ.infExp['infDiffcoef'] * cm2m
        nuFissionXS = univ.infExp['infNsf'] / cm2m
        try:
            Efiss = res.resdata['fisse'][step][0] * mev2j
        except:
            Efiss = res.resdata['fisse'][0] * mev2j
        powerXS = Efiss * univ.infExp['infFiss'] / cm2m
        scatterMatrixP0 = univ.infExp['infSp0'] / cm2m
        scatterMatrixP1 = univ.infExp['infSp1'] / cm2m
        scatterMatrixP2 = univ.infExp['infSp2'] / cm2m
        scatterMatrixP3 = univ.infExp['infSp3'] / cm2m
        scatterMatrixP4 = univ.infExp['infSp4'] / cm2m
        scatterMatrixP5 = univ.infExp['infSp5'] / cm2m

        scatterMatrixP0.shape = (ng, ng)
        scatterMatrixP1.shape = (ng, ng)
        scatterMatrixP2.shape = (ng, ng)
        scatterMatrixP3.shape = (ng, ng)
        scatterMatrixP4.shape = (ng, ng)
        scatterMatrixP5.shape = (ng, ng)

        totalXS = univ.infExp['infTot'] / cm2m
        removalXS = [tot - scatt[i] for i, tot, scatt in zip(range(ng), totalXS, scatterMatrixP0)]

        chiPrompt = univ.infExp['infChip']
        if (isEffective):
            chiDelayed = chiPrompt
            try:
                delayedFraction = res.resdata['adjIfpImpBetaEff'][step][::2]
            except:
                delayedFraction = res.resdata['adjIfpImpBetaEff'][::2]
        else:
            chiDelayed = univ.infExp['infChid']
            try:
                delayedFraction = res.resdata['fwdAnaBetaZero'][step][::2]
            except:
                delayedFraction = res.resdata['fwdAnaBetaZero'][::2]

        try:
            decayConstant = res.resdata['fwdAnaLambda'][step][::2]
        except:
            decayConstant = res.resdata['fwdAnaLambda'][::2]

        discFactor = [1] * ng

        univ0 = res.getUniv('0', burnup=burnup, index=step, timeDays=days)
        integralFlux = univ.infExp['infFlx'] / univ0.infExp['infFlx']

        # Fill the variables
        self.description = f"From Serpent universe '{univ.name}' ({burnup} MWd/kgU)"
        self.inverseVelocity = NonUniformList(inverseVelocity)
        self.diffusionCoefficient = NonUniformList(diffcoeff)
        self.nuFissionXS = NonUniformList(nuFissionXS)
        self.powerXS = NonUniformList(powerXS)
        self.removalXS = NonUniformList(removalXS)
        self.scatteringMatrixP0 = Matrix(scatterMatrixP0)
        self.scatteringMatrixP1 = Matrix(scatterMatrixP1)
        self.scatteringMatrixP2 = Matrix(scatterMatrixP2)
        self.scatteringMatrixP3 = Matrix(scatterMatrixP3)
        self.scatteringMatrixP4 = Matrix(scatterMatrixP4)
        self.scatteringMatrixP5 = Matrix(scatterMatrixP5)
        self.chiPrompt = NonUniformList(chiPrompt)
        self.chiDelayed = NonUniformList(chiDelayed)
        self.delayedFraction = NonUniformList(delayedFraction)
        self.decayConstant = NonUniformList(decayConstant)
        self.discFactor = NonUniformList(discFactor)
        self.integralFlux = NonUniformList(integralFlux)


    def read_from_openmc(
            self,
            outputFilename: str,
            domain,
            energyGroups,
            delayedGroups,
            rootDomain,
            openmcStatepointFile=None,
            scatteringMatrixNumber: int=6,
            isPrompt: bool=False
        ) -> None:
        """
        Extract from OpenMC's statepoint output file.

        Parameters
        ----------
        outputFilename : str
            OpenMC output file name (e.g `statepoint.1000.h5`)
        domain : openmc.Universe | openmc.Cell | openmc.Material
            Name of the OpenMC domain (universe, cell, material) to extract the
            MGXS
        energyGroups : openmc.mgxs.EnergyGroups
            Energy group structure used to create the MGXS
        delayedGroups
            List containing the neutron delay groups (e.g for 2 delayed groups
            [1, 2])
        rootDomain : openmc.Universe | openmc.Cell | openmc.Material
            Top level domain
        openmcStatepointFile : openmc.StatePoint
            OpenMC's StatePoint object containing the results
        scatteringMatrixNumber : int
            Number of scattering matrix order (default 6)
        isPrompt : bool
            Used for nuSigmaEff MGXS with `openmc.mgxs.FissionXS`
            (default False)
        """
        import openmc

        # Unit conversion from OpenMC to OpenFOAM
        ev2j = 1.602176487e-19
        cm2m = 0.01

        sp = openmcStatepointFile
        if (openmcStatepointFile is None):
            sp = openmc.StatePoint(outputFilename)

        ng = energyGroups.num_groups

        # Convert to MeV
        self.groupStructure = [e*1e-6 for e in energyGroups.group_edges]

        # Total XS
        totalXS = openmc.mgxs.TotalXS(
            domain=domain, energy_groups=energyGroups, name='Total'
        )

        # 1/V
        inverseVelocity = openmc.mgxs.InverseVelocity(
            domain=domain, energy_groups=energyGroups, name='IV'
        )

        # Diffision coefficient
        diffusionCoefficient = openmc.mgxs.DiffusionCoefficient(
            domain=domain, energy_groups=energyGroups, name='D'
        )

        # nuSigma_fission
        nuFissionXS = openmc.mgxs.FissionXS(
            domain=domain, energy_groups=energyGroups, nu=True,
            name="nuSigmaEff", prompt=isPrompt
        )

        # Sigma_power
        kappaFissionXS = openmc.mgxs.KappaFissionXS(
            domain=domain, energy_groups=energyGroups, name='sigmaPow'
        )

        # Scattering Matrix P0 to P5
        scatterMatrixXS = openmc.mgxs.ScatterMatrixXS(
            domain=domain, energy_groups=energyGroups,
            name='scatteringMatrixP', nu=True
        )
        scatterMatrixXS.legendre_order = scatteringMatrixNumber-1
        scatterMatrixXS.correction = None
        scatterMatrixXS.formulation = 'consistent'

        # Chi
        chiPrompt = openmc.mgxs.Chi(
            domain=domain, energy_groups=energyGroups, prompt=True,
            name='chiPrompt'
        )
        chiDelayed = openmc.mgxs.ChiDelayed(
            domain=domain, energy_groups=energyGroups,
            name='chiDelayed'
        )

        oneGroup = openmc.mgxs.EnergyGroups([0, 20e6])

        # Beta fraction
        beta = openmc.mgxs.Beta(
            domain=domain,
            energy_groups=oneGroup,
            delayed_groups=delayedGroups,
            name='Beta'
        )
        betaRoot = openmc.mgxs.Beta(
            domain=rootDomain,
            energy_groups=oneGroup,
            delayed_groups=delayedGroups,
            name='Beta'
        )

        # Decay rate constants
        Lambda = openmc.mgxs.DecayRate(
            domain=domain,
            energy_groups=oneGroup,
            delayed_groups=delayedGroups,
            name='lambda'
        )
        LambdaRoot = openmc.mgxs.DecayRate(
            domain=rootDomain,
            energy_groups=oneGroup,
            delayed_groups=delayedGroups,
            name='lambda'
        )

        totalXS.load_from_statepoint(sp)
        inverseVelocity.load_from_statepoint(sp)
        diffusionCoefficient.load_from_statepoint(sp)
        nuFissionXS.load_from_statepoint(sp)
        kappaFissionXS.load_from_statepoint(sp)
        scatterMatrixXS.load_from_statepoint(sp)
        chiPrompt.load_from_statepoint(sp)
        chiDelayed.load_from_statepoint(sp)
        beta.load_from_statepoint(sp)
        betaRoot.load_from_statepoint(sp)
        Lambda.load_from_statepoint(sp)
        LambdaRoot.load_from_statepoint(sp)
        integralFlux = sp.get_tally(name='integralFlux'+self.name).get_slice(scores=['flux'])
        fluxIntegralRoot = sp.get_tally(name='integralFluxmasterMultiGroup').get_slice(scores=['flux'])

        scatteringMatrix = lambda i: np.array([
            np.transpose(m)
            for m in np.transpose(
                scatterMatrixXS.get_xs()
            )
        ][i])


        P0_diag = [line[i] for i, line in enumerate(scatteringMatrix(0))]

        removalXS = np.array([
            inf_tot - P0_i for inf_tot, P0_i in zip(totalXS.get_xs(), P0_diag)
        ])

        # If all the lambda constant in the current MGXS object are at 0, use
        # the one provided by the master
        isCurrentLambdaZero = all(lambda_ == 0 for lambda_ in Lambda.get_xs())
        isCurrentBetaZero = all(beta_ == 0 for beta_ in beta.get_xs())

        decayConstant = LambdaRoot.get_xs() if isCurrentLambdaZero else Lambda.get_xs()
        delayedFraction = betaRoot.get_xs() if isCurrentBetaZero else beta.get_xs()

        discFactor = [1] * ng

        normalizedFlux = [
            flux[0][0]/fluxRoot[0][0] for flux, fluxRoot in zip(
                integralFlux.mean, fluxIntegralRoot.mean
            )
        ][::-1]

        self.description = f"From OpenMC universe '{domain.name}'"
        self.inverseVelocity = inverseVelocity.get_xs() / cm2m
        self.diffusionCoefficient = diffusionCoefficient.get_xs() * cm2m
        self.nuFissionXS = nuFissionXS.get_xs() / cm2m
        self.powerXS = kappaFissionXS.get_xs() * ev2j / cm2m
        self.removalXS = NonUniformList(removalXS / cm2m)
        self.scatteringMatrixP0 = Matrix(scatteringMatrix(0) / cm2m)
        self.scatteringMatrixP1 = Matrix(scatteringMatrix(1) / cm2m)
        self.scatteringMatrixP2 = Matrix(scatteringMatrix(2) / cm2m)
        self.scatteringMatrixP3 = Matrix(scatteringMatrix(3) / cm2m)
        self.scatteringMatrixP4 = Matrix(scatteringMatrix(4) / cm2m)
        self.scatteringMatrixP5 = Matrix(scatteringMatrix(5) / cm2m)
        self.chiPrompt = chiPrompt.get_xs()
        self.chiDelayed = chiDelayed.get_xs()[0]
        self.delayedFraction = NonUniformList([e[0] for e in delayedFraction])
        self.decayConstant = NonUniformList(decayConstant)
        self.discFactor = NonUniformList(discFactor)
        self.integralFlux = NonUniformList(normalizedFlux)


    def plot_xs_spectrum(
            self,
            ax,
            xsType: str,
            energyBins: list[float]=None,
            color: str=None,
            ls: str="-",
            label: str=None
        ):
        """
        Plot the XS in terms of energy discretization.

        Parameters
        ----------
        ax
            Matplotlib axis object to plot in.
        xsType : str
            Options: `sigmaRemoval`, `nuSigmaEff`, `sigmaPow`,
            `discFactor`, `chiPrompt`, `chiDelayed`, `IV`, `lambda`, `Beta`,
            `D`, `integralFlux`.
        energyBins : list[float]
            Energy discretization (default `None`). If `None`, use the group
            structure provided or find in Serpent/OpenMC. If group structure is
            `None`, use a linear binning of size one unit.
        color : str
            Line color (default `None`). If `None`, use the default colors of
            Matplotlib color wheel.
        ls : str
            Linestyle (default `"-"`).
        label : str
            Label of the line (default `None`). If `None`, use a default naming
            such as "{name} {xsType}".
        """
        check_type("xsType", xsType, str)
        check_value("xsType", xsType, _XS_TYPES)

        xsDict = {
            "sigmaRemoval": self.removalXS,
            "nuSigmaEff": self.nuFissionXS,
            "sigmaPow": self.powerXS,
            "discFactor": self.discFactor,
            "chiPrompt": self.chiPrompt,
            "chiDelayed": self.chiDelayed,
            "IV": self.inverseVelocity,
            "lambda": self.decayConstant,
            "Beta": self.delayedFraction,
            "D": self.diffusionCoefficient,
            "integralFlux": self.integralFlux,
        }

        xsToPlot = xsDict[xsType]

        if (self.groupStructure is not None):
            energyBins = self.groupStructure

        if (energyBins is None):
            energyBins = np.arange(-0.5, len(xsToPlot)+0.5, 1)

        xvalues, yvalues = [], []
        for ei, ef, xs in zip(energyBins[:-1], energyBins[1:], xsToPlot):
            xvalues.append(ei)
            xvalues.append(ef)
            yvalues.append(xs)
            yvalues.append(xs)

        if (label is None):
            label = f"{self.name} {xsType}"

        ax.plot(xvalues, yvalues, color=color, ls=ls, label=label)


    def plot_scattering_matrix(
            self,
            ax,
            xsType: str,
            energyBins: list[float]=None,
            isLogScale: bool=False,
            isSymLogScale: bool=False,
            zMin: float=None,
            zMax: float=None,
            cmap: str='RdBu_r',
            linthresh: float=None
        ):
        """
        Plot the scattering matrix XS in terms of energy discretization.

        Parameters
        ----------
        ax
            Matplotlib axis object to plot in.
        xsType : str
            Options: `scatteringMatrixP0`, `scatteringMatrixP1`,
            `scatteringMatrixP2`, `scatteringMatrixP3`, `scatteringMatrixP4`,
            `scatteringMatrixP5`.
        energyBins : list[float]
            Energy discretization (default `None`). If `None`, use the group
            structure provided or find in Serpent/OpenMC. If group structure is
            `None`, use a linear binning of size one unit.
        isLogScale : bool
            Log scale on the Z-axis (default `False`).
        isSymLogScale : bool
            Symmetric log scale on the Z-axis (default `False`).
        zMin : float
            Minimum Z, if `None`, compute the minimum value based on data
            (default `None`).
        zMax : float
            Minimum Z, if `None`, compute the maximum value based on data
            (default `None`).
        cmap : str
            Color map (default `RdBu_r`).
        linthresh : float
        """
        check_type("xsType", xsType, str)
        check_value("xsType", xsType, _XS_TYPES)

        import matplotlib.colors as colors

        xsDict = {
            "scatteringMatrixP0": self.scatteringMatrixP0,
            "scatteringMatrixP1": self.scatteringMatrixP1,
            "scatteringMatrixP2": self.scatteringMatrixP2,
            "scatteringMatrixP3": self.scatteringMatrixP3,
            "scatteringMatrixP4": self.scatteringMatrixP4,
            "scatteringMatrixP5": self.scatteringMatrixP5,
        }

        xsToPlot = np.array(xsDict[xsType])

        if (self.groupStructure is not None):
            energyBins = self.groupStructure

        if (energyBins is None):
            energyBins = np.arange(0, len(xsToPlot), 1)

        xRange, yRange = np.meshgrid(energyBins, energyBins)

        zMin = xsToPlot.flatten().min() if zMin is None else zMin
        zMax = xsToPlot.flatten().max() if zMax is None else zMax

        norm = None
        if (isLogScale):
            norm = colors.LogNorm(vmin=zMin, vmax=zMax)
        elif (isSymLogScale):
            norm = colors.SymLogNorm(vmin=zMin, vmax=zMax, base=10, linthresh=linthresh)
        else:
            norm = colors.Normalize(vmin=zMin, vmax=zMax)

        colormap = ax.pcolormesh(
            xRange, yRange, xsToPlot,
            cmap=cmap,
            norm=norm,
            shading='nearest'
        )
        ax.set_aspect('equal')
        ax.yaxis.set_inverted(True)

        return(colormap)


class NuclearDataState(OpenFOAMDict):
    """
    Nuclear data state object.
    The first Nuclear data state added to the NuclearData object must contain
    all the perturbation parameters if multiple state have been generated.

    Example ::

        refState = NuclearDataState(
            "reference",
            parameters={'Tfuel': 1000, 'rhoCool': 1000},
            zones=[...]
        )
        Tfuel1200K = NuclearDataState(
            "Tfuel1200K",
            parameters={'Tfuel': 1200}, # Assume by default that rhoCool is not perturb and uses rhoCool = 1000
            zones=[...]
        )
        TfuelAndRhoHot = NuclearDataState(
            "TfuelAndRhoHot",
            parameters={'Tfuel': 1200, 'rhoCool': 900},
            zones=[...]
        )

    Example of use with `read_from_openmc`::

        openmcState.read_from_openmc(
            outputFilename="../OpenMC/statepoint.500.h5",
            domains=[
                (openmc.Cell(name='Fuel'), 'fuel', 1),
                (openmc.Cell(name='Cladding'), 'cladding', 0),
                (openmc.Cell(name='Water'), 'water', 0),
            ],
            energyGroups=openmc.mgxs.EnergyGroups([0., 0.625, 20.0e6]),
            delayedGroups=list(range(1,7)),
            rootDomain=openmc.Universe(0, name='root universe')
        )

    Example of use with `read_from_serpent`::

        serpentState.read_from_serpent(
            outputFilename="../Serpent/main_res.m",
            universes=[
                ('uFuel', 'fuel', 1),
                ('uClad', 'cladding', 0),
                ('uWater', 'water', 0),
            ]
        )


    Parameters
    ----------
    name : str
        Name of the state
    parameters : dict
        Dictionary of physical parameters describing the state at which the
        nuclear data have been prepared.
        Example ::

            {'Tfuel': 300, 'rhoCool': 980}

    zones : list[NuclearDataZone]
        List of nuclear data zones
    """

    def __init__(
            self,
            name: str,
            parameters: dict[str, float]={},
            zones: list[NuclearDataZone]=[]
        ):
        super().__init__(name=name)

        self.zones: OpenFOAMList[NuclearDataZone] = OpenFOAMList(NuclearDataZone, "zones", items=zones)

        self.parameters = parameters


    def __repr__(self, depth = 0):
        for key, item in self.parameters.items():
            self.__setitem__(key, item)

        self.__setitem__("zones", self.zones)

        return super().__repr__(depth)


    def add_zone(self, zone: NuclearDataZone):
        check_type('zone', zone, NuclearDataZone)
        self.zones.append(zone)


    def add_parameters(self, name: str, value: int | float | Vector):
        """
        Add value at which the XS has been generated (e.g Tfuel = 1000).
        Can be called multiple times for multiple parameters (e.g Tfuel = 1000,
        rhoCool = 1000).
        """
        check_type('name', name, str)
        check_type('value', value, (int, float, Vector))
        self.parameters[name] = value


    def get_zone_by_name(self, zoneName: str) -> NuclearDataZone:
        for zone in self.zones:
            if (zoneName == zone.name):
                return(zone)

        return(None)


    def read_from_serpent(
            self,
            outputFilename: str,
            universes: list[tuple],
            step: int=0,
            burnup: float=None,
            days: float=None,
            isEffective: bool=True,
        ):
        """
        Extract the nuclear data from Serpent `.m` output file for multiple
        universes and a single burnup step.

        Example of use with `read_from_serpent`::

            serpentState.read_from_serpent(
                outputFilename="../Serpent/main_res.m",
                universes=[
                    ('uFuel', 'fuel', 1),
                    ('uClad', 'cladding', 0),
                    ('uWater', 'water', 0),
                ]
            )

        Parameters
        ----------
        outputFilename : str
            Serpent output file name
        universes : list[tuple]
            List of corresponding Serpent universe and OpenFOAM cellZone. Last
            number is the fuel fraction.
            Example ::

                [   # SSS2    OpenFOAM       fuelFraction
                    ('fuel', 'fuelAssembly', 0.55),
                    ('refl', 'reflector'   , 0),
                ]

        step : int
            Burnup step (default to 0)
        burnup : float
            Burnup [MWd/kgU] of the desired universe
        days : float
            Time [days] of the desired universe
        isEffective : bool
            Flag to extract the effective or analogue values (default to True)
        """
        import serpentTools

        # Share results for faster processing
        res = serpentTools.read(outputFilename)

        for universe, cellZone, fuelFraction in universes:
            newZone = NuclearDataZone(name=cellZone, fuelFraction=fuelFraction)
            newZone.read_from_serpent(
                outputFilename,
                universe=universe,
                burnup=burnup,
                step=step,
                days=days,
                isEffective=isEffective,
                serpentToolsResfile=res
            )
            self.add_zone(newZone)


    def read_from_openmc(
            self,
            outputFilename: str,
            domains: list[tuple],
            energyGroups,
            delayedGroups,
            rootDomain
        ):
        """
        Extract the nuclear data from OpenMC `statepoint.h5` output file for
        multiple domains.

        Example of use with `read_from_openmc`::

            openmcState.read_from_openmc(
                outputFilename="../OpenMC/statepoint.500.h5",
                domains=[
                    (openmc.Cell(name='Fuel'), 'fuel', 1),
                    (openmc.Cell(name='Cladding'), 'cladding', 0),
                    (openmc.Cell(name='Water'), 'water', 0),
                ],
                energyGroups=openmc.mgxs.EnergyGroups([0., 0.625, 20.0e6]),
                delayedGroups=list(range(1,7)),
                rootDomain=openmc.Universe(0, name='root universe')
            )

        Parameters
        ----------
        outputFilename : str
            OpenMC output file name
        domains : list[tuple]
            List of corresponding OpenMC domain and OpenFOAM cellZone. Last
            number is the fuel fraction. A domain can be `openmc.Universe`,
            `openmc.Cell`, or `openmc.Material`.

            Example ::

                [   # OpenMC                    OpenFOAM       fuelFraction
                    (openmc.Cell(name='fuel'), 'fuelAssembly', 0.55),
                    (openmc.Cell(name='refl'), 'reflector'   , 0),
                ]
        energyGroups : openmc.mgxs.EnergyGroups
            OpenMC energy groups used during the OpenMC simulation
        delayedGroups
            List containing the neutron delay groups (e.g for 2 delayed groups
            [1, 2])
        rootDomain : openmc.Universe | openmc.Cell | openmc.Material
            Top level domain
        """
        import openmc

        # Share results for faster processing
        with openmc.StatePoint(outputFilename) as sp:
            for domain, cellZone, fuelFraction in domains:
                newZone = NuclearDataZone(name=cellZone, fuelFraction=fuelFraction)
                newZone.read_from_openmc(
                    outputFilename,
                    domain=domain,
                    energyGroups=energyGroups,
                    delayedGroups=delayedGroups,
                    rootDomain=rootDomain,
                    openmcStatepointFile=sp,
                )
                self.add_zone(newZone)


    def plot_xs_spectrum(
            self,
            ax,
            zoneName: str,
            xsType: str,
            energyBins: list[float]=None,
            color: str=None,
            ls: str="-",
            label: str=None
        ):
        """
        Plot the XS in terms of energy discretization for all the zones
        contained in this state.

        Parameters
        ----------
        ax
            Matplotlib axis object to plot in.
        zoneName : str
            Name of the zone.
        xsType : str {"sigmaRemoval", "nuSigmaEff", "sigmaPow", "discFactor", "chiPrompt", "chiDelayed", "IV", "lambda", "Beta", "D", "integralFlux"}
            Nuclear data type.
        energyBins : list[float]
            Energy discretization (default `None`). If `None`, use the group
            structure provided or find in Serpent/OpenMC. If group structure is
            `None`, use a linear binning of size one unit.
        color : str
            Line color (default `None`). If `None`, use the default colors of
            Matplotlib color wheel.
        ls : str {"-", "--", "-.", ":"}
            Linestyle (default `"-"`)
        label : str
            Label of the line (default `None`). If `None`, use a default naming
            such as "{name} {xsType}".
        """
        for zone in self.zones:
            if (zone.name == zoneName):
                zone.plot_xs_spectrum(ax, xsType, energyBins, color, ls, label)
                break


class NuclearData(OpenFOAMFile):
    """
    Main Nuclear Data object. The nuclearData file contains the main parameter
    settings.

    Parameters
    ----------
    region : str
        Name of the region (default `""`)
    xsVariables : OpenFOAMDict
        Dictionary containing the XS variable laws (default `{}`)
    energyGroups : int
        Number of energy groups
    precGroups : int
        Number of precursor groups
    polyharmonicSplineMode : int
        Polyharmonic spline function mode (1. `r`; 2. `r^2 ln(r)`; 3. `r^3`;
        4. `r^4 ln(r)`) (default `1`).
    ScNo : float
        Schmidt number for diffusion of precursors (default `1`).
    adjustDiscFactors : bool
        Flag to apply the discontinuity factor adjustement (default `False`).
    useGivenDiscFactors : bool
        Flag to use homogeneous discontinuity factors provided in state/zone
        (default `False`)
    legendreMoments : int
        Number of Legendre moments, used in SN solver (default `None`).
    isLowMemory : bool
        Flag for low memory foot print, useful for SN solver (defaulf `None`).
    """
    def __init__(
            self,
            region: str="",
            xsVariables: OpenFOAMDict=None,
            energyGroups: int=None,
            precGroups: int=None,
            polyharmonicSplineMode: int=1,
            ScNo: float=1,
            axialOrientation: Vector=None,
            adjustDiscFactors: bool=False,
            useGivenDiscFactors: bool=False,
            legendreMoments: int=None,
            isLowMemory: bool=None,
        ):
        super().__init__("nuclearData", folder="constant", region=region)

        self.energyGroups = energyGroups
        self.precGroups = precGroups
        self.polyharmonicSplineMode = polyharmonicSplineMode
        self.ScNo = ScNo
        self.axialOrientation = axialOrientation
        self.adjustDiscFactors = adjustDiscFactors
        self.useGivenDiscFactors = useGivenDiscFactors
        self.legendreMoments = legendreMoments
        self.isLowMemory = isLowMemory

        self.xsVariables = xsVariables

        self.states = OpenFOAMList(NuclearDataState, "states")


    @property
    def energyGroups(self):
        return self._energyGroups

    @energyGroups.setter
    def energyGroups(self, energyGroups):
        check_type("energyGroups", energyGroups, (int, np.int64), none_ok=True)
        self._energyGroups = energyGroups

    @property
    def precGroups(self):
        return self._precGroups

    @precGroups.setter
    def precGroups(self, precGroups):
        check_type("precGroups", precGroups, (int, np.int64), none_ok=True)
        self._precGroups = precGroups

    @property
    def polyharmonicSplineMode(self):
        return self._polyharmonicSplineMode

    @polyharmonicSplineMode.setter
    def polyharmonicSplineMode(self, polyharmonicSplineMode):
        check_type("polyharmonicSplineMode", polyharmonicSplineMode, (int, np.int64), none_ok=True)
        self._polyharmonicSplineMode = polyharmonicSplineMode

    @property
    def ScNo(self):
        return self._ScNo

    @ScNo.setter
    def ScNo(self, ScNo):
        check_type("ScNo", ScNo, (float, int), none_ok=True)
        self._ScNo = ScNo

    @property
    def axialOrientation(self):
        return self._axialOrientation

    @axialOrientation.setter
    def axialOrientation(self, axialOrientation):
        check_type("axialOrientation", axialOrientation, (Vector, list, tuple), none_ok=True)
        if (isinstance(axialOrientation, (list, tuple))):
            self._axialOrientation = Vector(axialOrientation[0], axialOrientation[1], axialOrientation[2])
        else:
            self._axialOrientation = axialOrientation

    @property
    def adjustDiscFactors(self):
        return self._adjustDiscFactors

    @adjustDiscFactors.setter
    def adjustDiscFactors(self, adjustDiscFactors):
        check_type("adjustDiscFactors", adjustDiscFactors, bool)
        self._adjustDiscFactors = adjustDiscFactors

    @property
    def useGivenDiscFactors(self):
        return self._useGivenDiscFactors

    @useGivenDiscFactors.setter
    def useGivenDiscFactors(self, useGivenDiscFactors):
        check_type("useGivenDiscFactors", useGivenDiscFactors, bool)
        self._useGivenDiscFactors = useGivenDiscFactors

    @property
    def legendreMoments(self):
        return self._legendreMoments

    @legendreMoments.setter
    def legendreMoments(self, legendreMoments):
        check_type("legendreMoments", legendreMoments, int, none_ok=True)
        self._legendreMoments = legendreMoments

    @property
    def isLowMemory(self):
        return self._isLowMemory

    @isLowMemory.setter
    def isLowMemory(self, isLowMemory):
        check_type("isLowMemory", isLowMemory, bool, none_ok=True)
        self._isLowMemory = isLowMemory

    @property
    def xsVariables(self):
        return self._xsVariables

    @xsVariables.setter
    def xsVariables(self, xsVariables):
        check_type("xsVariables", xsVariables, OpenFOAMDict, none_ok=True)
        if (xsVariables is not None):
            self._xsVariables = xsVariables
        else:
            self._xsVariables = OpenFOAMDict({})

    # def add_external_file(self, srcpath: str) -> None:
    #     self.externalFiles.append(srcpath)


    def add_variable(self, name: str, law: str) -> None:
        """
        Add a variable law to XS perturbation parameter.

        Parameters
        ----------
        name : str
            Name of the variable
        law : str {"lin", "log", "sqrt"}
            Name of the law to be used by the interpolatation algorithm.
        """
        check_type("XS variable name", name, str)
        check_type("XS variable law", law, str)
        check_value("XS variable law", law, _XS_VARIABLE_LAW_TYPES)
        self.xsVariables[name] = law


    def add_state(self, state: NuclearDataState):
        check_type("state", state, NuclearDataState)
        self.states.append(state)


    def get_state_by_name(self, stateName: str) -> NuclearDataState:
        for state in self.states:
            if (state.name == stateName):
                return(state)

        return(None)


    def export_states_to_openfoam(self):
        text = "states\n(\n"
        for i, state in enumerate(self.states):
            # Replace the first state name to be the "reference" state
            if (i == 0):
                state.name = "reference"
            text += f"{tab}{state.__repr__(depth=1)}\n"

        text += ");\n"

        return(text)


    def export_for_spatial_solver_to_openfoam(self):
        ng = [zone.energyGroups for state in self.states for zone in state.zones]
        nd = [zone.precGroups for state in self.states for zone in state.zones]
        groupStructure = [
            zone.groupStructure
            for state in self.states
            for zone in state.zones
            if zone.groupStructure is not None
        ]

        text = ""

        if (len(ng) > 0):
            if (not all([e == ng[0] for e in ng])):
                ngErr = [e for e in ng if e != ng[0]][0]
                # idxErr = ng.index(ngErr)
                msg = "Inconsistency in the number of energy groups found in nuclear data zones.\n"
                msg += f"Found {ngErr} energy groups instead of {ng[0]}"
                raise ValueError(msg)

            if (not all([e == nd[0] for e in nd])):
                ndErr = [e for e in nd if e != nd[0]][0]
                msg = "Inconsistency in the number of precursor groups found in nuclear data zones.\n"
                msg += f"Found {ndErr} precursor groups instead of {nd[0]}"
                raise ValueError(msg)

            self.energyGroups = ng[0]
            self.precGroups = nd[0]

            if (len(groupStructure) > 0 and any([e != None for e in groupStructure[0]])):
                groupStructureStr = ', '.join([f'{e:g}' for e in groupStructure[0]])
                text += '// Group structure : [' + groupStructureStr + '] MeV\n'

            text += addParameter('energyGroups', self.energyGroups, isAddExtraLine=True)
            text += addParameter('precGroups', self.precGroups, isAddExtraLine=True)

            text += addParameter('polyharmonicSplineMode', self.polyharmonicSplineMode, isAddExtraLine=True)
            text += addParameter('ScNo', self.ScNo, isAddExtraLine=True)
            if (self.axialOrientation is not None):
                text += addParameter('axialOrientation', self.axialOrientation, isAddExtraLine=True)

            text += addParameter('adjustDiscFactors', self.adjustDiscFactors, isAddExtraLine=True, none_ok=False)
            text += addParameter('useGivenDiscFactors', self.useGivenDiscFactors, isAddExtraLine=True, none_ok=False)

            if (self.legendreMoments is not None):
                text += addParameter('legendreMoments', self.legendreMoments, isAddExtraLine=True, none_ok=False)
            if (self.isLowMemory is not None):
                text += addParameter('isLowMemory', self.isLowMemory, isAddExtraLine=True, none_ok=False)

        text += f"xsVariables{self.xsVariables!r}\n"

        text += self.export_states_to_openfoam()

        return(text)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self) -> str:
        text = ""
        text += self.export_for_spatial_solver_to_openfoam()
        return(text)


    def import_from_openfoam(self, path: str=None):
        """
        Import nuclear data from an already processed `nuclearData` file in the
        OpenFOAM format.
        """
        foamFile = foamlib.FoamFile(path if path is not None else self.path)

        self.energyGroups = foamFile['energyGroups']
        self.precGroups = foamFile['precGroups']
        if ('polyharmonicSplineMode' in foamFile):
            self.polyharmonicSplineMode = foamFile['polyharmonicSplineMode']
        if ('ScNo' in foamFile):
            self.ScNo = foamFile['ScNo']
        if ('axialOrientation' in foamFile):
            self.axialOrientation = foamFile['axialOrientation']
        if ('adjustDiscFactors' in foamFile):
            self.adjustDiscFactors = foamFile['adjustDiscFactors']
        if ('useGivenDiscFactors' in foamFile):
            self.useGivenDiscFactors = foamFile['useGivenDiscFactors']
        if ('legendreMoments' in foamFile):
            self.legendreMoments = foamFile['legendreMoments']
        if ('isLowMemory' in foamFile):
            self.isLowMemory = foamFile['isLowMemory']

        # Read variable law
        xsVariables = foamFile['xsVariables']
        for xsVariable in xsVariables:
            self.add_variable(name=xsVariable, law=xsVariables[xsVariable])

        # Read nuclear data states
        for state in foamFile["states"]:
            stateName, stateDict = state
            nuclearDataState = NuclearDataState(name=stateName)

            zones = stateDict['zones']
            for zone in zones:
                zoneName, values = zone

                nuclearDataZone = NuclearDataZone(
                    name=zoneName,
                    fuelFraction=values['fuelFraction'],
                    inverseVelocity=values['IV'],
                    diffusionCoefficient=values['D'],
                    nuFissionXS=values['nuSigmaEff'],
                    powerXS=values['sigmaPow'],
                    removalXS=values['sigmaRemoval'],
                    scatteringMatrixP0=values['scatteringMatrixP0'][1],
                    chiPrompt=values['chiPrompt'],
                    chiDelayed=values['chiDelayed'],
                    decayConstant=values['lambda'],
                    delayedFraction=values['Beta'],
                    discFactor=values['discFactor'],
                    integralFlux=values['integralFlux'],
                )
                if (f'scatteringMatrixP1' in values):
                    nuclearDataZone.scatteringMatrixP1 = values[f'scatteringMatrixP1'][1]
                if (f'scatteringMatrixP2' in values):
                    nuclearDataZone.scatteringMatrixP2 = values[f'scatteringMatrixP2'][1]
                if (f'scatteringMatrixP3' in values):
                    nuclearDataZone.scatteringMatrixP3 = values[f'scatteringMatrixP3'][1]
                if (f'scatteringMatrixP4' in values):
                    nuclearDataZone.scatteringMatrixP4 = values[f'scatteringMatrixP4'][1]
                if (f'scatteringMatrixP5' in values):
                    nuclearDataZone.scatteringMatrixP5 = values[f'scatteringMatrixP5'][1]

                nuclearDataState.add_zone(nuclearDataZone)

            self.add_state(nuclearDataState)



class PointKineticsData(NuclearData):
    """
    Special nuclear data object for point-kinetics simulation.
    """
    def __init__(
            self,
            region: str="",
            fastNeutrons: bool=None,
            promptGenerationTime: float=None,
            nuFission: float=None,
            energyPerFission: float=None,
            delayedFractions: list=[],
            decayConstants: list=[],
            feedbackCoeffDoppler: float=0,
            feedbackCoeffTFuel: float=0,
            feedbackCoeffTClad: float=0,
            feedbackCoeffTCool: float=0,
            feedbackCoeffRhoCool: float=0,
            feedbackCoeffTStruct: float=0,
            feedbackCoeffTStructMech: float=0,
            absoluteDrivelineExpansionCoeff: float=0,
            controlRodReactivityMap=None,
            externalReactivityTimeProfile: TimeProfile=None,
            boronReactivityTimeProfile: TimeProfile=None,
            fuelFeedbackZones: list=[],
            coolFeedbackZones: list=[],
            structFeedbackZones: list=[],
            drivelineFeedbackZones: list=[],
        ):
        super().__init__(region)

        self.fastNeutrons = fastNeutrons
        self.promptGenerationTime = promptGenerationTime
        self.nuFission = nuFission
        self.energyPerFission = energyPerFission

        self.delayedFractions = delayedFractions
        self.decayConstants = decayConstants

        self.feedbackCoeffDoppler = feedbackCoeffDoppler
        self.feedbackCoeffTFuel = feedbackCoeffTFuel
        self.feedbackCoeffTClad = feedbackCoeffTClad
        self.feedbackCoeffTCool = feedbackCoeffTCool
        self.feedbackCoeffRhoCool = feedbackCoeffRhoCool
        self.feedbackCoeffTStruct = feedbackCoeffTStruct
        self.feedbackCoeffTStructMech = feedbackCoeffTStructMech
        self.absoluteDrivelineExpansionCoeff = absoluteDrivelineExpansionCoeff

        self.controlRodReactivityMap: Table = controlRodReactivityMap

        self.externalReactivityTimeProfile: TimeProfile = externalReactivityTimeProfile
        self.boronReactivityTimeProfile: TimeProfile = boronReactivityTimeProfile

        self.fuelFeedbackZones = fuelFeedbackZones
        self.coolFeedbackZones = coolFeedbackZones
        self.structFeedbackZones = structFeedbackZones
        self.drivelineFeedbackZones = drivelineFeedbackZones

    @property
    def fastNeutrons(self):
        return self._fastNeutrons

    @fastNeutrons.setter
    def fastNeutrons(self, fastNeutrons) -> None:
        check_type("fastNeutrons", fastNeutrons, bool, none_ok=True)
        self._fastNeutrons = fastNeutrons

    @property
    def promptGenerationTime(self):
        return self._promptGenerationTime

    @promptGenerationTime.setter
    def promptGenerationTime(self, promptGenerationTime) -> None:
        check_type("promptGenerationTime", promptGenerationTime, (float, int), none_ok=True)
        self._promptGenerationTime = promptGenerationTime

    @property
    def nuFission(self):
        return self._nuFission

    @nuFission.setter
    def nuFission(self, nuFission) -> None:
        check_type("nuFission", nuFission, (float, int), none_ok=True)
        self._nuFission = nuFission

    @property
    def energyPerFission(self):
        return self._energyPerFission

    @energyPerFission.setter
    def energyPerFission(self, energyPerFission) -> None:
        check_type("energyPerFission", energyPerFission, (float, int), none_ok=True)
        self._energyPerFission = energyPerFission

    @property
    def delayedFractions(self):
        return self._delayedFractions

    @delayedFractions.setter
    def delayedFractions(self, delayedFractions) -> None:
        check_type("delayedFractions", delayedFractions, (list, np.ndarray), none_ok=True)
        self._delayedFractions = delayedFractions

    @property
    def decayConstants(self):
        return self._decayConstants

    @decayConstants.setter
    def decayConstants(self, decayConstants) -> None:
        check_type("decayConstants", decayConstants, (list, np.ndarray), none_ok=True)
        self._decayConstants = decayConstants

    @property
    def feedbackCoeffDoppler(self):
        return self._feedbackCoeffDoppler

    @feedbackCoeffDoppler.setter
    def feedbackCoeffDoppler(self, feedbackCoeffDoppler) -> None:
        check_type("feedbackCoeffDoppler", feedbackCoeffDoppler, (float, int), none_ok=True)
        self._feedbackCoeffDoppler = feedbackCoeffDoppler

    @property
    def feedbackCoeffTFuel(self):
        return self._feedbackCoeffTFuel

    @feedbackCoeffTFuel.setter
    def feedbackCoeffTFuel(self, feedbackCoeffTFuel) -> None:
        check_type("feedbackCoeffTFuel", feedbackCoeffTFuel, (float, int), none_ok=True)
        self._feedbackCoeffTFuel = feedbackCoeffTFuel

    @property
    def feedbackCoeffTClad(self):
        return self._feedbackCoeffTClad

    @feedbackCoeffTClad.setter
    def feedbackCoeffTClad(self, feedbackCoeffTClad) -> None:
        check_type("feedbackCoeffTClad", feedbackCoeffTClad, (float, int), none_ok=True)
        self._feedbackCoeffTClad = feedbackCoeffTClad

    @property
    def feedbackCoeffTCool(self):
        return self._feedbackCoeffTCool

    @feedbackCoeffTCool.setter
    def feedbackCoeffTCool(self, feedbackCoeffTCool) -> None:
        check_type("feedbackCoeffTCool", feedbackCoeffTCool, (float, int), none_ok=True)
        self._feedbackCoeffTCool = feedbackCoeffTCool

    @property
    def feedbackCoeffRhoCool(self):
        return self._feedbackCoeffRhoCool

    @feedbackCoeffRhoCool.setter
    def feedbackCoeffRhoCool(self, feedbackCoeffRhoCool) -> None:
        check_type("feedbackCoeffRhoCool", feedbackCoeffRhoCool, (float, int), none_ok=True)
        self._feedbackCoeffRhoCool = feedbackCoeffRhoCool

    @property
    def feedbackCoeffTStruct(self):
        return self._feedbackCoeffTStruct

    @feedbackCoeffTStruct.setter
    def feedbackCoeffTStruct(self, feedbackCoeffTStruct) -> None:
        check_type("feedbackCoeffTStruct", feedbackCoeffTStruct, (float, int), none_ok=True)
        self._feedbackCoeffTStruct = feedbackCoeffTStruct

    @property
    def feedbackCoeffTStructMech(self):
        return self._feedbackCoeffTStructMech

    @feedbackCoeffTStructMech.setter
    def feedbackCoeffTStructMech(self, feedbackCoeffTStructMech) -> None:
        check_type("feedbackCoeffTStructMech", feedbackCoeffTStructMech, (float, int), none_ok=True)
        self._feedbackCoeffTStructMech = feedbackCoeffTStructMech

    @property
    def absoluteDrivelineExpansionCoeff(self):
        return self._absoluteDrivelineExpansionCoeff

    @absoluteDrivelineExpansionCoeff.setter
    def absoluteDrivelineExpansionCoeff(self, absoluteDrivelineExpansionCoeff) -> None:
        check_type("absoluteDrivelineExpansionCoeff", absoluteDrivelineExpansionCoeff, (float, int), none_ok=True)
        self._absoluteDrivelineExpansionCoeff = absoluteDrivelineExpansionCoeff

    @property
    def controlRodReactivityMap(self):
        return self._controlRodReactivityMap

    @controlRodReactivityMap.setter
    def controlRodReactivityMap(self, controlRodReactivityMap) -> None:
        check_type("controlRodReactivityMap", controlRodReactivityMap, (list, np.ndarray, Table), none_ok=True)
        if (isinstance(controlRodReactivityMap, Table)):
            self._controlRodReactivityMap = controlRodReactivityMap
            self._controlRodReactivityMap.type = 'controlRodReactivityMap'
        elif (controlRodReactivityMap is None):
            self._controlRodReactivityMap = Table([], type='controlRodReactivityMap')
        else:
            self._controlRodReactivityMap = Table(controlRodReactivityMap, type='controlRodReactivityMap')

    @property
    def fuelFeedbackZones(self):
        return self._fuelFeedbackZones

    @fuelFeedbackZones.setter
    def fuelFeedbackZones(self, fuelFeedbackZones) -> None:
        check_type("fuelFeedbackZones", fuelFeedbackZones, (list, np.ndarray), none_ok=True)
        self._fuelFeedbackZones = fuelFeedbackZones

    @property
    def coolFeedbackZones(self):
        return self._coolFeedbackZones

    @coolFeedbackZones.setter
    def coolFeedbackZones(self, coolFeedbackZones) -> None:
        check_type("coolFeedbackZones", coolFeedbackZones, (list, np.ndarray), none_ok=True)
        self._coolFeedbackZones = coolFeedbackZones

    @property
    def structFeedbackZones(self):
        return self._structFeedbackZones

    @structFeedbackZones.setter
    def structFeedbackZones(self, structFeedbackZones) -> None:
        check_type("structFeedbackZones", structFeedbackZones, (list, np.ndarray), none_ok=True)
        self._structFeedbackZones = structFeedbackZones

    @property
    def drivelineFeedbackZones(self):
        return self._drivelineFeedbackZones

    @drivelineFeedbackZones.setter
    def drivelineFeedbackZones(self, drivelineFeedbackZones) -> None:
        check_type("drivelineFeedbackZones", drivelineFeedbackZones, (list, np.ndarray), none_ok=True)
        self._drivelineFeedbackZones = drivelineFeedbackZones


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter("fastNeutrons", self.fastNeutrons, isAddExtraLine=True)
        text += addParameter("promptGenerationTime", self.promptGenerationTime, isAddExtraLine=True)

        if (self.nuFission is not None):
            text += addParameter("nuFission", self.nuFission, isAddExtraLine=True)
        if (self.energyPerFission is not None):
            text += addParameter("energyPerFission", self.energyPerFission, isAddExtraLine=True)

        text += addParameter("Beta", List(self.delayedFractions), isAddExtraLine=True)
        text += addParameter("lambda", List(self.decayConstants), isAddExtraLine=True)

        text += addParameter("feedbackCoeffDoppler", self.feedbackCoeffDoppler, isAddExtraLine=True)
        text += addParameter("feedbackCoeffTFuel", self.feedbackCoeffTFuel, isAddExtraLine=True)
        text += addParameter("feedbackCoeffTClad", self.feedbackCoeffTClad, isAddExtraLine=True)
        text += addParameter("feedbackCoeffTCool", self.feedbackCoeffTCool, isAddExtraLine=True)
        text += addParameter("feedbackCoeffRhoCool", self.feedbackCoeffRhoCool, isAddExtraLine=True)
        text += addParameter("feedbackCoeffTStruct", self.feedbackCoeffTStruct, isAddExtraLine=True)
        text += addParameter("feedbackCoeffTStructMech", self.feedbackCoeffTStructMech, isAddExtraLine=True)
        text += addParameter("absoluteDrivelineExpansionCoeff", self.absoluteDrivelineExpansionCoeff, isAddExtraLine=True)

        text += f"{self.controlRodReactivityMap!r};\n\n"

        if (self.externalReactivityTimeProfile is not None):
            text += f"externalReactivityTimeProfile{self.externalReactivityTimeProfile!r}\n"

        if (self.boronReactivityTimeProfile is not None):
            text += f"boronReactivityTimeProfile{self.boronReactivityTimeProfile!r}\n"

        text += addParameter("fuelFeedbackZones", List(self.fuelFeedbackZones), isAddExtraLine=True)
        text += addParameter("coolFeedbackZones", List(self.coolFeedbackZones), isAddExtraLine=True)
        text += addParameter("structFeedbackZones", List(self.structFeedbackZones), isAddExtraLine=True)
        text += addParameter("drivelineFeedbackZones", List(self.drivelineFeedbackZones), isAddExtraLine=True)

        text += self.export_for_spatial_solver_to_openfoam()

        return(text)


class ControlRodMovement(OpenFOAMDict):
    def __init__(
            self,
            cellZoneName: str,
            startTime: float,
            endTime: float,
            speed: float,
            initialDistanceFromMeshCR: float,
            followerName: str
        ):
        super().__init__(name=cellZoneName)

        self.startTime = startTime
        self.endTime = endTime
        self.speed = speed
        self.initialDistanceFromMeshCR = initialDistanceFromMeshCR
        self.followerName = followerName

    @property
    def startTime(self):
        return self._startTime

    @startTime.setter
    def startTime(self, startTime) -> None:
        check_type("startTime", startTime, (float, int))
        self._startTime = startTime
        self.__setitem__("startTime", startTime)

    @property
    def endTime(self):
        return self._endTime

    @endTime.setter
    def endTime(self, endTime) -> None:
        check_type("endTime", endTime, (float, int))
        self._endTime = endTime
        self.__setitem__("endTime", endTime)

    @property
    def speed(self):
        return self._speed

    @speed.setter
    def speed(self, speed) -> None:
        check_type("speed", speed, (float, int))
        self._speed = speed
        self.__setitem__("speed", speed)

    @property
    def initialDistanceFromMeshCR(self):
        return self._initialDistanceFromMeshCR

    @initialDistanceFromMeshCR.setter
    def initialDistanceFromMeshCR(self, initialDistanceFromMeshCR) -> None:
        check_type("initialDistanceFromMeshCR", initialDistanceFromMeshCR, (float, int))
        self._initialDistanceFromMeshCR = initialDistanceFromMeshCR
        self.__setitem__("initialDistanceFromMeshCR", initialDistanceFromMeshCR)

    @property
    def followerName(self):
        return self._followerName

    @followerName.setter
    def followerName(self, followerName) -> None:
        check_type("followerName", followerName, str)
        self._followerName = followerName
        self.__setitem__("followerName", followerName)



class ControlRodMove(OpenFOAMFile):
    def __init__(
            self,
            region = ""
        ):
        super().__init__("CRmove", "constant", region)

        self.zones = OpenFOAMList(ControlRodMovement, "zones")


    def add_control_rod_movement(
            self,
            cellZoneName: str,
            startTime: float,
            endTime: float,
            speed: float,
            initialDistanceFromMeshCR: float,
            followerName: str
        ):
        check_type("cellZoneName", cellZoneName, str)

        self.zones.append(ControlRodMovement(
            cellZoneName=cellZoneName,
            startTime=startTime,
            endTime=endTime,
            speed=speed,
            initialDistanceFromMeshCR=initialDistanceFromMeshCR,
            followerName=followerName
        ))

    @property
    def is_empty(self) -> bool:
        return(self.zones.is_empty)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += f"{self.zones!r}\n"

        return(text)
