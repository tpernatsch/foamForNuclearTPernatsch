"""
Script to automatically generate the externalSourceFlux{i} files for GeN-Foam
external source simulations.
Each externalSourceFlux file can process an external PTrack file to generate the
external source. Each file is generated for one energy bin.

Author: Thomas Guilbaud, 2022/11/11
"""

#------------------------------------ Main -------------------------------------

# energies=[0.0, 8.0, 20.0]
# energies=[1.9640330000E+01,
#     1.0000000000E+01,
#     6.0653070000E+00,
#     3.6787940000E+00,
#     2.2313020000E+00,
#     1.3533530000E+00,
#     8.2085000000E-01,
#     4.9787070000E-01,
#     3.0197380000E-01,
#     1.8315640000E-01,
#     1.1109000000E-01,
#     6.7379470000E-02,
#     4.0867710000E-02,
#     2.4787520000E-02,
#     1.5034390000E-02,
#     9.1188200000E-03,
#     5.5308440000E-03,
#     3.3546260000E-03,
#     2.0346840000E-03,
#     1.2340980000E-03,
#     7.4851830000E-04,
#     4.5399930000E-04,
#     3.0432480000E-04,
#     1.4862540000E-04,
#     9.1660880000E-05,
#     6.7904050000E-05,
#     4.0169000000E-05,
#     2.2603290000E-05,
#     1.3709590000E-05,
#     8.3152870000E-06,
#     4.0000000000E-06,
#     5.4000000000E-07,
#     1.0000000000E-07,
#     1.0000100000E-11
# ]
# energies.reverse()
# energies = [0.0, 1e-4, 0.03, 0.8, 20.0]
# energies = [1e-9, 0.01, 0.2, 1.8, 20.0]
energies = [20.0, 1.8, 0.2, 0.01, 1e-9]

beamIntensity = 0.003366917 # 5.06216e-3 # 1.51316e-3 # 5.252e-4 # 3.67e-3 # 4.39e-3 # A

nu_s = 21.154 # neutrons/proton, spallation yield

# -------------------------------- End of main --------------------------------

for i, lowEnergy, highEnergy in zip(range(len(energies)-1), energies[1:], energies[:-1]):
    file = open(f"externalSourceFlux{i}", 'w')
    file.write("""/*--------------------------------*- C++ -*----------------------------------*\\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \\  /  |/ /  ______  / /_     / __ \\ / __ `/  / __ `__ \\    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \\____/   \\___/ /_/ |_/          /_/       \\____/ \\__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      defaultExternalSourceFlux;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 -3 -1 0 0 0 0 ];

// internalField uniform 0;

internalField   #codeStream
{
    code
    #{
        const IOdictionary& d = static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(d.db());

        // Set field
        scalarField sourceFluxField(mesh.nCells(), 0.0); // Default to 0
        scalarField volume(mesh.V());

        // Read file
        std::ifstream infile(\"../storedData/0_ptrack.out\");
        // Set boundary values
        const double energyMin = """+str(lowEnergy)+", energyMax = "+str(highEnergy)+""";
        const double xMin = -0.7, xMax = 0.7;
        const double yMin = -0.7, yMax = 0.7;
        const double zMin = -1.0, zMax = 1.52;

        Info<< "Energy range = [" << energyMin << "; " << energyMax << "] MeV"
            << endl;

        // Set temp values during reading
        double x, y, z, dx, dy, dz, energyI, energyO, wgt, timePart;
        double count = 0;

        // Loop over the lines of the PTrack file
        while (infile >> x >> y >> z >> dx >> dy >> dz >> energyI >> energyO >> wgt >> timePart)
        {
            x *= 1e-2; y *= 1e-2; z = z*1e-2 + 0.8;
            if (xMin <= x && x <= xMax && yMin <= y && y <= yMax && zMin <= z && z <= zMax)
            {
                if (energyMin <= energyO && energyO < energyMax)
                {
                    label cellNumber = mesh.findCell(point(x, y, z));
                    if (cellNumber >= 0)
                    {
                        sourceFluxField[cellNumber] += 1.0 / volume[cellNumber];
                    }
                }
                count += 1.0;
            }
        }

        // Normalise the neutron source with the proton beam current and neutron yield
        const double intensity = """+str(beamIntensity)+""", nu_s = """+str(nu_s)+""", eV = 1.602e-19;
        const double scale = nu_s * intensity / eV; // nb neutrons/s
        Info<< \"Scaling factor = \" << scale << \" neutrons/s\" << endl;
        forAll(sourceFluxField, sourceFluxFieldI)
        {
            sourceFluxField[sourceFluxFieldI] *= scale / count;
        }

        // Write the new internalField
        sourceFluxField.writeEntry(\"\", os);
    #};

    //! Optional:
    codeInclude
    #{
        #include <fstream>
        #include \"fvCFD.H\"
        #include \"meshSearch.H\"
    #};

    //! Optional:
    codeOptions
    #{
        -I$(LIB_SRC)/finiteVolume/lnInclude \\
        -I$(LIB_SRC)/meshTools/lnInclude
    #};
};

boundaryField
{
    ".*"
    {
        type            fixedValue;
        value           uniform 0;
    }
}

// ************************************************************************* //
""")
    file.close()
