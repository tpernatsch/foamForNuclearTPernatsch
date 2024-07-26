"""
Title: Multi-Group XS per cell for OpenMC. Allow to convert the generated MGXS
from OpenMC into GeN-Foam format.

Author: Thomas Guilbaud, 2022/07/17

Sources:
    https://docs.openmc.org/en/v0.13.0/pythonapi/mgxs.html
"""
#==============================================================================*
#                                   Imports
#==============================================================================*

import numpy as np
import openmc

#==============================================================================*
#                              Usefull Functions
#==============================================================================*

def printVector(file, name: str, vector: list, scale: float=1) -> None:
    """
        Function to print a vector inline for GeN-Foam nuclearData file.

        name: Name of the vector to print.
        vector: 1D-array to print.
        scale: Scaling factor applied to all the values in the vector.
    """
    file.write("\t\t{:18} nonuniform List<scalar> {} ({});\n".format(
        name,
        len(vector),
        " ".join(["{:.6e}".format(e*scale) for e in vector])
    ))

def printMatrix(file, name: str, matrix: list, scale: float=1) -> None:
    """
        Function to print a matrix inline for GeN-Foam nuclearData file.

        name: Name of the matrix to print.
        matrix: 2D-array to print.
        scale: Scaling factor applied to all the values in the matrix.
    """
    file.write("\t\t{0} {1} {1} (\n".format(
        name, len(matrix)
    ))
    for line in matrix:
        temp = "\t\t\t("
        temp += " ".join(["{:.6e}".format(e*scale) for e in line])
        temp += ")\n"
        file.write(temp)
    file.write("\t\t);\n")


#==============================================================================*
#                                  Main Class
#==============================================================================*

class MultiGroupXS(object):
    """docstring for MultiGroupXS."""

    def __init__(
        self,
        domain,
        energyGroups: openmc.mgxs.EnergyGroups,
        delayedGroups: list,
        fuelFraction: float=0,
        genfoamName: str="0",
        isPrompt: bool=False
    ) -> None:
        """
            Constructor of the Multi-Group Cross-Section object.

            domain: Cell, Universe, or Material.
            energyGroups: List containing the boundaries of the energy groups
                (e.g [0, 0.625, 20e6]).
            delayedGroups: List containing the neutron delay groups (e.g for 2
                delayed groups [1, 2]).
            fuelFraction: Fraction of fuel inside the domain.
            genfoamName: Name for the GeN-Foam nuclearData.zones object.
            isPrompt: Used for nuSigmaEff MGXS with openmc.mgxs.FissionXS.
        """
        self.domain = domain
        self.energyGroups = energyGroups
        oneGroupeBoundaries = np.array(
            [energyGroups.group_edges[0], energyGroups.group_edges[-1]]
        )
        self.oneGroup = openmc.mgxs.EnergyGroups(group_edges=oneGroupeBoundaries)
        self.delayedGroups = delayedGroups
        self.isPrompt = isPrompt

        self.sp = None

        self.genfoamName = genfoamName
        self.fuelFraction = fuelFraction
        self.scatteringMatrixNumber = 6
        self.totalXS = None
        self.inverseVelocity = None
        self.diffusionCoefficient = None
        self.nuFissionXS = None
        self.kappaFissionXS = None
        self.scatterMatrixXS = None
        self.chiPrompt = None
        self.chiDelayed = None
        self.beta = None
        self.Lambda = None
        self.integralFlux = None
        self.promptGenerationTime = 0

        self.GenerateTallies()


    def AddToTallies(self, tallies: openmc.Tallies) -> None:
        """
            Method to add all the MGXS specific tallies to the master tally
            object used to generate the XML.

            tallies: openmc.Tallies object to store the MGXS tallies.
        """
        tallies += self.totalXS.tallies.values()
        tallies += self.inverseVelocity.tallies.values()
        tallies += self.diffusionCoefficient.tallies.values()
        tallies += self.nuFissionXS.tallies.values()
        tallies += self.kappaFissionXS.tallies.values()
        tallies += self.scatterMatrixXS.tallies.values()
        tallies += self.chiPrompt.tallies.values()
        tallies += self.chiDelayed.tallies.values()
        tallies += self.beta.tallies.values()
        tallies += self.Lambda.tallies.values()
        tallies.append(self.integralFlux)


    def GenerateTallies(self) -> None:
        """
            Generate and allocate the MGXS specific tallies for GeN-Foam
            nuclearData files.
        """

        # Total XS
        self.totalXS = openmc.mgxs.TotalXS(
            domain=self.domain, energy_groups=self.energyGroups, name='Total'
        )

        # 1/V
        self.inverseVelocity = openmc.mgxs.InverseVelocity(
            domain=self.domain, energy_groups=self.energyGroups, name='IV'
        )

        # Diffision coefficient
        self.diffusionCoefficient = openmc.mgxs.DiffusionCoefficient(
            domain=self.domain, energy_groups=self.energyGroups, name='D'
        )

        # nuSigma_fission
        self.nuFissionXS = openmc.mgxs.FissionXS(
            domain=self.domain, energy_groups=self.energyGroups, nu=True,
            name="nuSigmaEff", prompt=self.isPrompt
        )

        # Sigma_power
        self.kappaFissionXS = openmc.mgxs.KappaFissionXS(
            domain=self.domain, energy_groups=self.energyGroups, name='sigmaPow'
        )

        # Scattering Matrix P0 to P5
        self.scatterMatrixXS = openmc.mgxs.ScatterMatrixXS(
            domain=self.domain, energy_groups=self.energyGroups,
            name='scatteringMatrixP'
        )
        self.scatterMatrixXS.legendre_order = self.scatteringMatrixNumber-1
        self.scatterMatrixXS.correction = None
        self.scatterMatrixXS.formulation = 'consistent'

        # Chi
        self.chiPrompt = openmc.mgxs.Chi(
            domain=self.domain, energy_groups=self.energyGroups, prompt=True,
            name='chiPrompt'
        )
        self.chiDelayed = openmc.mgxs.ChiDelayed(
            domain=self.domain, energy_groups=self.energyGroups,
            name='chiDelayed'
        )

        # Beta fraction
        self.beta = openmc.mgxs.Beta(
            domain=self.domain,
            energy_groups=self.oneGroup,
            delayed_groups=self.delayedGroups,
            name='Beta'
        )

        # Decay rate constants
        self.Lambda = openmc.mgxs.DecayRate(
            domain=self.domain,
            energy_groups=self.oneGroup,
            delayed_groups=self.delayedGroups,
            name='lambda'
        )

        # Integral flux over the domain
        self.integralFlux = openmc.Tally(name='integralFlux'+self.genfoamName)
        energyFilter = openmc.EnergyFilter(self.energyGroups.group_edges)
        self.integralFlux.filters = [energyFilter]
        if (type(self.domain) == openmc.Cell):
            self.integralFlux.filters.append(openmc.CellFilter(self.domain))
        elif (type(self.domain) == openmc.Universe):
            self.integralFlux.filters.append(openmc.UniverseFilter(self.domain))
        elif (type(self.domain) == openmc.Material):
            self.integralFlux.filters.append(openmc.MaterialFilter(self.domain))
        self.integralFlux.scores = ['flux', 'inverse-velocity']


    def GetPromptGenerationLifetime(self, isFluxWeighted: bool=True):
        """
            Return the physical prompt generation lifetime using the integrated
            inverse velocity tally.
        """
        # https://openmc.discourse.group/t/kinetic-parameter-calculation/596
        # There is currently (Version 0.13.0) no way to compute the EFFECTIVE
        # generation lifetime in OpenMC. The 1/v gives the flux-weighted prompt
        # removal lifetime.

        if (isFluxWeighted):
            self.promptGenerationTime = self.sp.get_tally(
                name='integralFlux'+self.genfoamName).get_slice(
                    scores=['inverse-velocity'])
            return(
                self.promptGenerationTime.mean[0][0][0],
                self.promptGenerationTime.std_dev[0][0][0]
            )
        
        # Other method using 1/(v * nu * Sigma_f)
        IV = self.inverseVelocity.get_xs()[0]
        nuFiss = self.nuFissionXS.get_xs()[0]
        promptGenerationTime = IV / nuFiss
        relErrIV = self.inverseVelocity.get_xs(value='rel_err')[0]
        relErrNuFiss = self.nuFissionXS.get_xs(value='rel_err')[0]
        promptGenerationTimeError = promptGenerationTime * np.sqrt(relErrIV**2 + relErrNuFiss**2)

        return(promptGenerationTime, promptGenerationTimeError)


    def LoadFromStatepoint(self, sp: openmc.StatePoint) -> None:
        """
            Load from StatePoint file the MGXS.

            sp: StatePoint containing the results of the simulation.
        """
        self.sp = sp
        self.totalXS.load_from_statepoint(sp)
        self.inverseVelocity.load_from_statepoint(sp)
        self.diffusionCoefficient.load_from_statepoint(sp)
        self.nuFissionXS.load_from_statepoint(sp)
        self.kappaFissionXS.load_from_statepoint(sp)
        self.scatterMatrixXS.load_from_statepoint(sp)
        self.chiPrompt.load_from_statepoint(sp)
        self.chiDelayed.load_from_statepoint(sp)
        self.beta.load_from_statepoint(sp)
        self.Lambda.load_from_statepoint(sp)
        self.integralFlux = sp.get_tally(
            name='integralFlux'+self.genfoamName
        ).get_slice(scores=['flux'])


    def PrintForGeNFoam(self, file, fluxIntegralMaster, Lambda: list=[]) -> None:
        """
            Print the MGXS results in the GeN-Foam format for the domain.

            file: File where the results needs to be written.
            Lambda: List of the decay constants for all the precursors groups.
        """

        # Usefull constants
        nGroups = self.energyGroups.num_groups
        cm2m = 1e-2
        ev2j = 1.602176487e-19

        # Write for nuclearData file
        file.write("\t{}  // OpenMC: {}\n".format(self.genfoamName, self.domain.name))
        file.write("\t{\n")
        file.write("\t\tfuelFraction       {:.6e};\n".format(self.fuelFraction))

        printVector(file, self.inverseVelocity.name, self.inverseVelocity.get_xs(), 1/cm2m)

        printVector(file, self.diffusionCoefficient.name, self.diffusionCoefficient.get_xs(), cm2m)

        printVector(file, self.nuFissionXS.name, self.nuFissionXS.get_xs(), 1/cm2m)

        printVector(file, self.kappaFissionXS.name, self.kappaFissionXS.get_xs(), ev2j/cm2m)

        scatteringMatrix = lambda i: [
            np.transpose(m) for m in np.transpose(self.scatterMatrixXS.get_xs())][i]
        for i in range(self.scatteringMatrixNumber):
            printMatrix(
                file,
                self.scatterMatrixXS.name+str(i),
                scatteringMatrix(i),
                1/cm2m
            )

        P0_diag = [line[i] for i, line in enumerate(scatteringMatrix(0))]
        sigmaDiapp = list(map(
            lambda inf_tot, P0_i: inf_tot - P0_i, self.totalXS.get_xs(), P0_diag
        ))
        printVector(file, "sigmaDisapp", sigmaDiapp, 1/cm2m)

        printVector(file, self.chiPrompt.name, self.chiPrompt.get_xs())

        printVector(file, self.chiDelayed.name, self.chiDelayed.get_xs()[0])

        printVector(file, self.beta.name, [e[0] for e in self.beta.get_xs()])

        # If all the lambda constant in the current MGXS object are at 0, use
        # the one provided by the master
        isCurrentLambdaZero = all(lambda_ == 0 for lambda_ in self.Lambda.get_xs())
        printVector(
            file,
            self.Lambda.name,
            Lambda.get_xs() if isCurrentLambdaZero else self.Lambda.get_xs()
        )

        printVector(file, "discFactor", [1]*nGroups)

        normalizedFlux = [
            flux[0][0]/fluxMaster[0][0] for flux, fluxMaster in zip(
                self.integralFlux.mean, fluxIntegralMaster.mean
            )
        ]
        printVector(file, "integralFlux", normalizedFlux)

        file.write("\t}\n\n")

#==============================================================================*
