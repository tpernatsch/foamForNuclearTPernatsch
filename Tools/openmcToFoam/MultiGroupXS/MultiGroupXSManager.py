"""
Title: Multi-Group XS manager for OpenMC. Allow to convert the generated MGXS
from OpenMC into GeN-Foam format.

Author: Thomas Guilbaud, 2022/07/17
"""

#==============================================================================*
#                                   Imports
#==============================================================================*

from MultiGroupXS import *


#==============================================================================*
#                                  Main Class
#==============================================================================*

class MultiGroupXSManager(object):
    """docstring for MultiGroupXSManager."""

    def __init__(
        self,
        domain,
        energyGroups: openmc.mgxs.EnergyGroups,
        delayedGroups: list
    ) -> None:
        """
            Constructor of the Multi-Group Cross-Section manager object.

            domain: Master Cell, Universe, or Material.
            energyGroups: List containing the boundaries of the energy groups
                (e.g [0, 0.625, 20e6]).
            delayedGroups: List containing the neutron delay groups (e.g for 2
                delayed groups [1, 2]).
        """
        self.mgxsList = []
        self.energyGroups = energyGroups
        self.delayedGroups = delayedGroups
        self.masterMGXSOneGroup = MultiGroupXS(
            domain,
            openmc.mgxs.EnergyGroups([0., 20.0e6]), # equivalent to one-group
            delayedGroups,
            genfoamName='masterOneGroup',
            isPrompt=True
        )
        self.masterMGXS = MultiGroupXS(
            domain,
            energyGroups,
            delayedGroups,
            genfoamName='masterMultiGroup',
            isPrompt=True
        )

        self.AddDomain(self.masterMGXSOneGroup)
        self.AddDomain(self.masterMGXS)


    def AddDomain(self, MultiGroupXS: MultiGroupXS) -> None:
        """
            Add a MGXS domain to be computed for GeN-Foam.

            MultiGroupXS: MGXS to be computed.
        """
        self.mgxsList.append(MultiGroupXS)


    def AddToTallies(self, tallies: openmc.Tallies) -> None:
        """
            Method to add all the MGXS specific tallies to the master tally
            object used to generate the XML.

            tallies: openmc.Tallies object to store the MGXS tallies.
        """
        for mgxsCell in self.mgxsList:
            mgxsCell.AddToTallies(tallies)


    def GenerateTallies(self) -> None:
        """
            Generate and allocate all the MGXS specific tallies for GeN-Foam
            nuclearData files.
        """
        for mgxsCell in self.mgxsList:
            mgxsCell.GenerateTallies()


    def LoadFromStatepoint(self, sp: openmc.StatePoint) -> None:
        """
            Load from StatePoint file the MGXS for all MGXS object stored.

            sp: StatePoint containing the results of the simulation.
        """
        for mgxsCell in self.mgxsList:
            mgxsCell.LoadFromStatepoint(sp)


    def PrintForGeNFoam(self, filename: str="nuclearData") -> None:
        """
            Print all the informations of the MGXS data for GeN-Foam.
            Print various standard and mandatory inputs for the user.

            filename: Name of the GeN-Foam resulting file (e.g nuclearData).
        """
        file = open(filename, "w")

        file.write(r"/*--------------------------------*- C++ -*----------------------------------*\ ")
        file.write(r"|       ______          _   __           ______                               |")
        file.write(r"|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |")
        file.write(r"|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |")
        file.write(r"|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |")
        file.write(r"|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |")
        file.write(r"|    Copyright (C) 2015 - 2022 EPFL                                           |")
        file.write(r"\*---------------------------------------------------------------------------*/")
        file.write("FoamFile\n")
        file.write("{\n")
        file.write("\tversion     2.0;\n")
        file.write("\tformat      ascii;\n")
        file.write("\tclass       dictionary;\n")
        file.write('\tlocation    "constant/neutroRegion";\n')
        file.write("\tobject      {};\n".format(filename))
        file.write("}\n")
        file.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")

        file.write("energyGroups                {}; // [{}] eV\n\n".format(
            self.energyGroups.num_groups, 
            ", ".join([f"{val:.6e}" for val in self.energyGroups.group_edges])
        ))
        file.write("precGroups                  {};\n\n".format(len(self.delayedGroups)))

        Lambda = self.masterMGXSOneGroup.Lambda

        if (filename == "nuclearData"):
            file.write("//- Point Kinetics data\n")
            file.write("fastNeutrons                false;\n\n")

            # Prompt Generation Lifetime
            promptGenerationTime, promptGenerationTimeError = self.masterMGXSOneGroup.GetPromptGenerationLifetime()
            file.write("promptGenerationTime        {:.6e}; // [s] ({:.6e})\n\n".format(promptGenerationTime, promptGenerationTimeError))

            # Beta
            file.write("Beta\n(\n")
            beta = self.masterMGXSOneGroup.beta
            betaTot = 0
            for mean, std_dev in zip(beta.get_xs(), beta.get_xs(value='std_dev')):
                file.write("\t{:.6e}  // ({:.6e})\n".format(mean[0], std_dev[0]))
                betaTot += mean[0]
            file.write(");\n")
            file.write("// Beta tot = {:.3f} pcm\n\n".format(betaTot*1e5))

            # Lambda
            file.write("lambda\n(\n")
            for mean, std_dev in zip(Lambda.get_xs(), Lambda.get_xs(value='std_dev')):
                file.write("\t{:.6e}  // ({:.6e})\n".format(mean, std_dev))
            file.write(");\n\n")

            # Feedback coefficients
            file.write("feedbackCoeffFastDoppler    0;\n\n")
            file.write("//- Representative of fuel axial expansion\n")
            file.write("feedbackCoeffTFuel          0; // [1/K]\n\n")
            file.write("//- Representative of in-assembly structure density change\n")
            file.write("//  (i.e. cladding AND wrapper wire)\n")
            file.write("feedbackCoeffTClad          0; // [1/K]\n\n")
            file.write("feedbackCoeffTCool          0; // [1/K]\n\n")
            file.write("feedbackCoeffRhoCool        0; // [1/(kg/m3)]\n\n")
            file.write("//- Representative of structure radial expansion\n")
            file.write("feedbackCoeffTStruct        0; // [1/K]\n\n")
            file.write("absoluteDrivelineExpansionCoeff 0; // [m/K]\n\n")
            file.write("//- CR insertion/driveline expansion (in m) vs reactivity\n")
            file.write("controlRodReactivityMap ();\n\n\n")

            # Zones
            file.write("fuelFeedbackZones ();\n\n")
            file.write("coolantFeedbackZones ();\n\n")
            file.write("structuresFeedbackZones ();\n\n")
            file.write("drivelineFeedbackZones ();\n\n\n")

        # Diffusion data
        file.write("//- Diffusion data\n")
        file.write("zones\n(\n")
        # Print all the zones except the first two, which are the master
        # MultiGroupXS for one-group and multi-group created in the constructor
        # above
        for mgxsCell in self.mgxsList[2:]:
            mgxsCell.PrintForGeNFoam(
                file,
                fluxIntegralMaster=self.masterMGXS.integralFlux,
                Lambda=Lambda
            )
        file.write(");\n\n")

        file.write("// ************************************************************************* //\n")

        file.close()

#==============================================================================*
