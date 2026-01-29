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

#include "hardening.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "Tuple2.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hardening, 0);
    addToRunTimeSelectionTable
    (
        yieldStressModel, 
        hardening, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hardening::hardening
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    yieldStressModel(mesh, lawDict)
{
#ifdef OPENFOAMFOUNDATION
    stressPlasticStrainSeries_.set(new Function1s::Table<scalar>
        (
            "plasticStrainVsYieldStress", lawDict
        )
    );
#elif OPENFOAMESI
    stressPlasticStrainSeries_.reset(new Table
        (
            "plasticStrainVsYieldStress", lawDict
        )
    );
#endif    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hardening::~hardening()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hardening::correctYieldStress
(
    const labelList& addr
)
{  

    const volScalarField& epsilonPEq_(mesh_.lookupObject<volScalarField>("epsilonPEq"));
    const scalarField& epsilonPEqI = epsilonPEq_.internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();

    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];
    
        sigmaY_[cellI] = stressPlasticStrainSeries_->value(epsilonPEqI[cellI]);
        scalar sigmaYOld = stressPlasticStrainSeries_->value(epsilonPEqOldI[cellI]);

        // Define linear plastic modulus
        Hp_[cellI] = (sigmaY_[cellI] - sigmaYOld)/max(
            epsilonPEqI[cellI] - epsilonPEqOldI[cellI], SMALL);

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                scalarField& sigmaYP(sigmaY_.boundaryFieldRef()[patchID]);
                scalarField& HpP(Hp_.boundaryFieldRef()[patchID]);
                
                const scalarField& epsilonPEqP
                (
                    epsilonPEq_.boundaryField()[patchID]
                );
                
                const scalarField& epsilonPEqOldP
                (
                    epsilonPEq_.oldTime().boundaryField()[patchID]
                );

                sigmaYP[faceID] = stressPlasticStrainSeries_->value(epsilonPEqP[faceID]);
                scalar sigmaYPOld = stressPlasticStrainSeries_->value(epsilonPEqOldP[faceID]);
                HpP[faceID] = (sigmaYP[faceID] - sigmaYPOld)/max(
                    epsilonPEqP[faceID] - epsilonPEqOldP[faceID], SMALL);
            }
        }
    }
}
// ************************************************************************* //