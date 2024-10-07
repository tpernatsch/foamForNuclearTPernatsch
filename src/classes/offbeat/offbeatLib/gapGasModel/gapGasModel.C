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

#include "gapGasModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gapGasModel, 0);
    defineRunTimeSelectionTable(gapGasModel, dictionary);

    addToRunTimeSelectionTable
    (
        gapGasModel, 
        gapGasModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::gapGasModel::correctMassFractions()
{
    scalar Ytot = sum(Y_);
    
    if (Ytot < SMALL)
    {
        FatalErrorIn
        (
            "Foam::gapFRAPCON::correctMassFractions"
        ) << "Sum of mass fractions is zero or negative"
          << exit(FatalError);
    }
    
    forAll(Y_, i)
    {
        Y_[i] /= Ytot;
    }
}


void Foam::gapGasModel::correctMolFractions()
{    
    scalar molTot(0.0);
    
    forAll(M_, i)
    {
        M_[i] = Y_[i]*gasM_/speciesW_[i];
        molTot += M_[i];
    }
    
    if (molTot < SMALL)
    {
        FatalErrorIn
        (
            "Foam::gapFRAPCON::correctMolFractions"
        ) << "Sum of molar fractions is zero or negative"
          << exit(FatalError);
    }
    
    M_ = M_/molTot;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gapGasModel::gapGasModel
(
    const fvMesh& mesh,
    const materials& mat,
    const volScalarField& T,
    const volVectorField& DorDD,
    const fissionGasRelease& fgr,
    const dictionary& gapGasOptDict
)
:
    regIOobject
    (
        IOobject
        (
            "gapGas",
            mesh.time().timeName()/"uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh),
    mat_(mat),
    T_(T),
    DorDD_(DorDD),
    fgr_(fgr),
    gapGasOptDict_(gapGasOptDict),
    species_(),
    speciesW_(),
    gapV_(0),
    gapT_(290),
    gasP_(0),
    gasM_(0),
    gasM0_(0),
    Y_(),
    Y0_(),
    M_(),
    M0_()
{
}

Foam::autoPtr<Foam::gapGasModel>
Foam::gapGasModel::New
(   
    const fvMesh& mesh,
    const materials& mat,
    const volScalarField& T,
    const volVectorField& DorDD,
    const fissionGasRelease& fgr,
    const dictionary& solverDict
)
{
    // Initialize default type for gapGasModel class
    word type("gapFRAPCON");

    dictionary gapGasOptDict
    (
        solverDict.subOrEmptyDict("gapGasOptions")
    );

    solverDict.lookup("gapGas") >> type;
    Info << "Selecting gapGasModel " << type << endl;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "gapGasModel::New\
            (const fvMesh&, const materials&, const volScalarField&, const volVectorField&, const fissionGasRelease&, const dictionary&)"
        )
            << "Unknown gapGasModel type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting gapGasModel type "
            << type << endl;
    }

    return autoPtr<gapGasModel>
    (
        cstrIter()
        (       
            mesh,  
            mat, 
            T, 
            DorDD,
            fgr,
            gapGasOptDict
        )
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gapGasModel::~gapGasModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::gapGasModel::updateVariables()
{
    gasM0_ = gasM_;
    Y0_ = Y_;
    M0_ = M_;
}


void Foam::gapGasModel::correctMass()
{  
    // Obtain list of gas components included in the fgr model
    const wordList& gasNames(fgr_.gasComponents());

    // Obtain mol of gas released as calculated by fgr model
    const scalarList& deltaGasMol(fgr_.gasMols());

    // Rescale mass fraction list using old time values
    Y_ = Y0_*gasM0_;

    // Rescale current gas mass using old time value
    gasM_ = gasM0_;

    // Loop over all gas species 
    forAll(gasNames, gasNameI)
    {
            const word gasName = gasNames[gasNameI];
            bool found(false);
            
            forAll(species_, specieI)
            {
                    if(species_[specieI] == gasName)
                    {
                            Y_[specieI] += 
                            deltaGasMol[gasNameI]*speciesW_[specieI]/1000;
                            
                            gasM_ += 
                            deltaGasMol[gasNameI]*speciesW_[specieI]/1000;
                            
                            found = true;
                            break;
                    }
            }

            if(!found)
            {
                FatalErrorIn("gapGasModel::correctMass")
                << "Gas component "
                << gasNames[gasNameI] << nl 
                << " not found in gapFRAPCON" << abort(FatalError);
            }
    }  

    correctMassFractions();
    correctMolFractions();     
}

// ************************************************************************* //
