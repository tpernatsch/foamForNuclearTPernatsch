/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fuelMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fuelMaterial, 0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::dictionary Foam::fuelMaterial::readIsotopesComposition() 
{

    if(not materialModelDict_.found("isotopes"))
    {
        // Lookup enrichment 
        scalar e(readScalar(materialModelDict_.lookup("enrichment")));

        // Create isotopes dictionary
        dictionary d("isotopes");

        // Create the sub-dictionary for Uranium
        dictionary uraniumDict;
        uraniumDict.add("ratioOverMetal", 1.0);
        uraniumDict.add("massNumbers", List<label>{235, 238});
        uraniumDict.add("weightFractions", List<scalar>{e, 1 - e});

        // Add the sub-dictionary under the key "U"
        d.add("U", uraniumDict);

        return d;
    }
    else
    {
        return materialModelDict_.subDict("isotopes");
    }
}

void Foam::fuelMaterial::checkIsotopesComposition() 
const
{
    // Get list of metals constituent of fuel
    const wordList fuelMetals(isotopesDict_.toc());
    
    // Init sum of concentrations of all the constituents (U, Pu, Np, etc.)
    scalar sumWeightConcentrations(0);

    // Loop over all fuel constituents
    Info<< tab << " Reading mass fractions of heavy metals in fuel..." << nl;
    forAll(fuelMetals, i)
    {
        const word metalName(fuelMetals[i]);
        const dictionary& metalDict(isotopesDict_.subDict(metalName));
        sumWeightConcentrations += readScalar(metalDict.lookup("ratioOverMetal"));
        
        // Lookup list of isotopes' mass numbers
        const List<scalar> massNumbers
        (
            metalDict.lookup("massNumbers")
        );

        // Lookup list of isotopes' weight fractions
        const List<scalar> weightFractions
        (
            metalDict.lookup("weightFractions")
        );

        if (massNumbers.size() != weightFractions.size())
        {
            FatalErrorInFunction
                << "Inconsistent isotope specification for U: "
                << "massNumbers size = " << massNumbers.size()
                << ", weightFractions size = " << weightFractions.size()
                << exit(FatalError);
        }
        
        // Init sum of isotopic concentrations for this element
        scalar sumIsotopes(0);

        Info<< tab << " " << metalName << " - " 
        << readScalar(metalDict.lookup("ratioOverMetal"))*100 << " %wt : " << nl;

        // Loop over all isotopes for this element
        forAll(massNumbers, j)
        {
            const scalar A_j = massNumbers[j];
            const scalar w_j = weightFractions[j];

            Info<< tab << tab 
                << " isotope MM " << A_j << " g/mol : " << w_j*100 << " %wt"
                << endl;

            sumIsotopes += w_j;
        }

        // Throw error if sum of element's isotopes != 1 
        if
        ( 
            sumIsotopes > 1.001 
            or 
            sumIsotopes < 0.999
        )
        {
            FatalErrorIn("Foam::fuelMaterial::fuelMaterial()")
            << "The sum of the concentrations for the element " << metalName  
            << " is " << sumIsotopes << ", but it must be equal to 1." << nl 
            << nl << exit(FatalError);
        }

    }
    
    // Throw error if sum of constituents' concentrations != 1 
    if
    ( 
        sumWeightConcentrations > 1.001 
        or 
        sumWeightConcentrations < 0.999
    )
    {
        FatalErrorIn("Foam::fuelMaterial::fuelMaterial()")
        << "The sum of the concentrations of the provided nuclides is " <<
        sumWeightConcentrations << ", but it must be equal to 1." << nl << nl
        << exit(FatalError);
    }
}

void Foam::fuelMaterial::initFuelIsotopesField() 
{
    // Get list of metals constituent of fuel
    const wordList list_metalsName(isotopesDict_.toc());

    // 
    scalarList list_MM_HM(list_metalsName.size(), 0.0);

    // Define avogadro number
    const scalar nAvo(6.022e23);

    // Define molar mass oxygen in g/mol;
    const scalar MMO(16);

    // Compute density of this material
    const scalar TD(readScalar(materialModelDict_.lookup("theoreticalDensity")));
    const scalar rho(densityFrac_ * TD);

    // Compute molar mass of heavy metals
    scalar MM_HM(0);
    
    // Loop over all fuel metal constituents
    forAll(list_metalsName, i)
    {
        const word metalName(list_metalsName[i]);
        const dictionary& metalDict(isotopesDict_.subDict(metalName));
        const scalar metalWeightFrac(readScalar(metalDict.lookup("ratioOverMetal")));
        
        // Lookup list of isotopes' mass numbers
        const List<scalar> massNumbers
        (
            metalDict.lookup("massNumbers")
        );

        // Lookup list of isotopes' weight fractions
        const List<scalar> weightFractions
        (
            metalDict.lookup("weightFractions")
        );

        if (massNumbers.size() != weightFractions.size())
        {
            FatalErrorInFunction
                << "Inconsistent isotope specification for U: "
                << "massNumbers size = " << massNumbers.size()
                << ", weightFractions size = " << weightFractions.size()
                << exit(FatalError);
        }

        forAll(massNumbers, j)
        {
            const scalar A_j = massNumbers[j];
            const scalar w_j = weightFractions[j];

            list_MM_HM[i] += A_j * w_j;
        }

        // Add to MM of fuel compound
        MM_HM += metalWeightFrac * list_MM_HM[i];
    }

    // Add oxygen MM to MM fuel (only metals have been added for the moment)
    const scalar MM_compound = MM_HM + OM_ * MMO;

    // Compute number of moles per cm3 of the compound (HM)O2 - this is equal to number of moles of HM
    const scalar n_HM = (rho * 1000 / 1e6) / MM_compound;

    // Compute number density of metallic nuclides (1/cm3) 
    const scalar N_HM = n_HM * nAvo;

    // Compute number density of oxigen nuclides (1/cm3) - this is "OM times" the number density of HM
    const scalar N_O = OM_ * N_HM;

    // Now create fields for each of this nuclides
    forAll(list_metalsName, i)
    {
        const word metalName(list_metalsName[i]);
        const dictionary& metalDict(isotopesDict_.subDict(metalName));

        // weight frac of this metal over all heavy metals
        const scalar metalWeightFrac(readScalar(metalDict.lookup("ratioOverMetal")));
        
        // Lookup list of isotopes' mass numbers
        const List<scalar> massNumbers
        (
            metalDict.lookup("massNumbers")
        );

        // Lookup list of isotopes' weight fractions
        const List<scalar> weightFractions
        (
            metalDict.lookup("weightFractions")
        );

        if (massNumbers.size() != weightFractions.size())
        {
            FatalErrorInFunction
                << "Inconsistent isotope specification for U: "
                << "massNumbers size = " << massNumbers.size()
                << ", weightFractions size = " << weightFractions.size()
                << exit(FatalError);
        }

        // Add number density of this element (sum of all isotopes) to registry
        volScalarField& N_element = 
        createOrLookup<scalar>
        (
            mesh_, 
            "N_" + metalName, 
            dimless/dimVolume, 
            0.0, 
            "zeroGradient"
        );

        // Covert mass fraction to molar fraction for this element
        const scalar x_el = metalWeightFrac * MM_HM / list_MM_HM[i];

        // init number density of this element
        const scalar N_el = x_el * N_HM;

        // Assign number densities of nuclides per cell (1/cm3)
        forAll(addr_, k)
        {
            const label cellI(addr_[k]);
            N_element[cellI] = N_el;
        }

        forAll(massNumbers, j)
        {
            const scalar A_j = massNumbers[j];
            const scalar w_j = weightFractions[j];

            // Convert mass fraction to molar fraction for this isotope
            const scalar x_iso = w_j * list_MM_HM[i] / A_j;

            // Add this nuclide to registry
            volScalarField& N_isotope = 
            createOrLookup<scalar>
            (
                mesh_, 
                "N_" + metalName + Foam::name(A_j), 
                dimless/dimVolume, 
                0.0, 
                "zeroGradient"
            );

            // Assign number densities of nuclides per cell (1/cm3)
            forAll(addr_, k)
            {
                const label cellI(addr_[k]);
                N_isotope[cellI] = x_iso * N_el;
            }
        }
    }

    // Add also field for oxygen
    volScalarField& Noxygen = 
    createOrLookup<scalar>
    (
        mesh_, 
        "N_O16", 
        dimless/dimVolume, 
        0.0, 
        "zeroGradient"
    );
    // Assign number densities of nuclides per cell (1/cm3)
    forAll(addr_, k)
    {
        const label cellI(addr_[k]);
        Noxygen[cellI] = N_O;
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fuelMaterial::fuelMaterial
(
    const fvMesh& mesh, 
    const dictionary& materialModelDict,
    const labelList& addr
)
:
    materialModel(mesh, materialModelDict, addr),
    isotopesDict_(readIsotopesComposition()),
    densification_(),
    swelling_(),
    relocation_(),
    failure_(),
    poreVelocity_(),
    rGrain_(readScalar(materialModelDict_.lookup("rGrain"))),
    OM_(materialModelDict_.lookupOrDefault("oxygenMetalRatio", 2.0)),
    oxygenMetalRatio_
    (
        createOrLookup<scalar>
        (
            mesh, 
            "oxygenMetalRatio", 
            dimless, 
            0.0,
            "zeroGradient"
        )
    ),
    densityFrac_(readScalar(materialModelDict_.lookup("densityFraction"))),
    porosity_
    (
        createOrLookup<scalar>
        (
            mesh, 
            "porosity", 
            dimless, 
            0.0,
            "zeroGradient"
        )
    ),
    dishFraction_(readScalar(materialModelDict_.lookup("dishFraction", 0.0)))
{
    // Initialize porosity and oxygen to metal ratio
    forAll(addr_, i)
    {
        const label cellI = addr_[i];   

        porosity_[cellI] = 1 - densityFrac_;
        oxygenMetalRatio_[cellI] = OM_;
    }
 
    porosity_.correctBoundaryConditions();

    // Check consistent input of isotopes compositions
    checkIsotopesComposition();

    // Initialize Isotopes Fields [1/cm3]
    initFuelIsotopesField();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fuelMaterial::~fuelMaterial()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
