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

#include "relocationFrapcon.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relocationFrapcon, 0);
    addToRunTimeSelectionTable(relocationModel, relocationFrapcon, dictionary);

    const char* relocationFrapcon::group_ = 
        ("behavioralModels::relocation::UO2Frapcon");

    defineParameter(relocationFrapcon, F_epsilonRelocation_, 
        "F_epsilonRelocation", (dimless), 1.0);
        
    defineParameter(relocationFrapcon, delta_epsilonRelocation_, 
        "delta_epsilonRelocation", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relocationFrapcon::relocationFrapcon
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    relocationModel(mesh, dict, baseDict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    heatSourceName_(dict.lookupOrDefault<word>("heatSourceName", "Q" )),
    gapWidthName_(dict.lookupOrDefault<word>("gapWidthName", "gapWidth" )),
    outerPatchID_(0),
    gapWidth_(nullptr),
    mapper_(nullptr), 
    GapCold_(0),
    DiamCold_(0),
    recFrac_(dict.lookupOrDefault<scalar>("recoveryFraction", 0.0)),
    relaxRecovery_(dict.lookupOrDefault<scalar>("relaxRecovery", 1.0)),    
    modifiedRelocationModel_(dict.lookupOrDefault<Switch>("modifiedRelocationModel", true))
{
    dict.lookup("GapCold") >> GapCold_;
    dict.lookup("DiamCold") >> DiamCold_;

    word outerPatch(dict.lookup("outerPatch"));
    outerPatchID_ = mesh_.boundaryMesh().findPatchID(outerPatch);
    
    if (outerPatchID_ == -1)
    {
        FatalIOErrorInFunction(dict.lookup("outerPatch"))
            << "Specified patch " << outerPatch << " not found."
            << abort(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relocationFrapcon::~relocationFrapcon()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relocationFrapcon::correct
(
    const scalarField& T, 
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F_relocation = 
    isInUserParameters(this->group(), "F_epsilonRelocation") ?
    this->F_epsilonRelocation() : 
    relocationModel::F_epsilonRelocation();

    const dimensionedScalar& F_relocationRecovery = 
    isInUserParameters(this->group(), "F_epsilonRelocationRecovery") ?
    this->F_epsilonRelocationRecovery() : 
    relocationModel::F_epsilonRelocationRecovery();

    const dimensionedScalar& delta_relocation = 
    isInUserParameters(this->group(), "delta_epsilonRelocation") ?
    this->delta_epsilonRelocation() : 
    relocationModel::delta_epsilonRelocation();

    const dimensionedScalar& delta_relocationRecovery = 
    isInUserParameters(this->group(), "delta_epsilonRelocationRecovery") ?
    this->delta_epsilonRelocationRecovery() : 
    relocationModel::delta_epsilonRelocationRecovery();
        
    epsilonRelocation_.storeOldTimes();
    epsilonRecoveredRelocation_.storeOldTimes();
    
    // Correct only when Q, Bu and gWidth are properly initialized
    if
    (
        mesh_.foundObject<volScalarField>(heatSourceName_) and
        mesh_.foundObject<volScalarField>(burnupName_) and
        mesh_.foundObject<volScalarField>(gapWidthName_)and
        mesh_.foundObject<sliceMapper>("sliceMapper")
    )
    { 
        if( gapWidth_ == nullptr )
        {
            gapWidth_ = &mesh_.lookupObject<volScalarField>(gapWidthName_);
        }

        if(mapper_ == nullptr)
        {
            mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
        }

        scalar nSlices(mapper_->nSlices());
        if(nSlices<=0)
        {
            FatalErrorInFunction()
            << "Number of slices cannot be less than 1." << nl
            << "Check that a sliceMapper different than \"none\" is used "
            << "or check the slice definition."
            << abort(FatalError);
        }

        // Take const references to heat source and burnup fields
        const volScalarField& Q = mesh_.lookupObject<volScalarField>(heatSourceName_);
        const volScalarField& Bu = mesh_.lookupObject<volScalarField>(burnupName_);

        // Calculate the slice-average values of Bu and Q 
        const scalarField& Q_1D(mapper_->sliceAverage(Q));
        const scalarField& Bu_1D(mapper_->sliceAverage(Bu));

        scalarField gWidth_1D(nSlices, 0.0);
        
        const scalarField& Vi = mesh_.V();

        const fvPatchScalarField& gapWidthP = 
        gapWidth_->boundaryField()[outerPatchID_];

        // NOTE: the following loop is not optimal as each material model
        // will calculate averages for every slice in the fuel.
        // (only a problem for simulation with multiple fuel materials)
        const PtrList<labelList>& sliceAddrList(mapper_->sliceAddrList());
        const List<bool>& isFuelList(mapper_->isFuel());

        forAll(sliceAddrList, i)
        {
            if(isFuelList[i])
            {
                const labelList& sliceAddr = sliceAddrList[i];

                scalarField sliceV(Vi, sliceAddr);

                List<cell> sliceCells(mesh_.cells(), sliceAddr);
                scalarField sliceGapWidth(sliceAddr.size(), 0.0);
                
                scalar minGapWidth(GREAT);

                forAll(sliceCells, j)
                {
                    const cell c(sliceCells[j]);

                    forAll(c, faceI)
                    {   
                        const label patchID =
                        mesh_.boundaryMesh().whichPatch(c[faceI]);

                        if (patchID == outerPatchID_)
                        {
                            const label faceID =
                            mesh_.boundaryMesh()[outerPatchID_].whichFace(c[faceI]);

                            const scalar faceSf =
                            gapWidthP.patch().magSf()[faceID];

                            if (gapWidthP[faceID] < minGapWidth)
                            {
                                minGapWidth = gapWidthP[faceID];
                            }
                        }  
                    }

                }

                gWidth_1D[i] = minGapWidth;
            }
        }

        reduce(gWidth_1D, minOp<scalarField>());

        const scalarField& sliceIDs(mapper_->sliceID().internalField());

        // Reference to relocation, old time relocation and recovered relocation
        scalarField& rel = epsilonRelocation_.ref();
        const scalarField& relOld = epsilonRelocation_.oldTime().internalField();
        scalarField& relRecover = epsilonRecoveredRelocation_.ref();
        const scalarField& relRecoverOld = epsilonRecoveredRelocation_.oldTime().internalField();

        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];
            const label sliceID = sliceIDs[cellI];

            // Calculate relocation. Relocation cannot decrease!
            scalar newRelocation = setRelocation
                                    (
                                        Q_1D[sliceID],
                                        GapCold_,
                                        DiamCold_,
                                        Bu_1D[sliceID]
                                    );

            // Perturb epsilon Relocation for sensitivity analysis
            newRelocation = newRelocation * F_relocation.value() 
                        + delta_relocation.value();

            scalar relocation = max(relOld[cellI], newRelocation);

            // if(rel[cellI]==relOld[cellI])
            // {
            //     relRecover[cellI] = 0;
            // }
            
            // Surface penetration
            scalar penetration = 
                ((-gWidth_1D[sliceID] + relRecover[cellI]*DiamCold_/2.0)  > SMALL) ? 
                (-gWidth_1D[sliceID] + relRecover[cellI]*DiamCold_/2.0) : 0;

            // // Relocation cannot cause stress!
            // if(penetration >= relocation*DiamCold_/2 )
            // {
            //               _______________________________ _ _ _ _ _ _ _ _       
            //                |    /|\           /|\          |   /|\  
            //                |     |penetration  |R_fuel     |    |relocation*R_fuel
            //                |     |             |        _ _|_ _\|/_ _ _ _ _
            //        ________|_____|_____________|___________|_______________
            //       /|\      |                   |           |
            //        |       |                   |           |
            //        |Ri_clad|                   |           |
            //        |       |                   |           |
                

            //     // Relocation is recovered as much as possible to prevent penetration
            //     rel[cellI] = relocation*(1 - recFrac_);
            //     relRecover[cellI] = relocation*(recFrac_);
            // }
            // else
            {
                /*          _______________________________ _ _ _ _ _ _ _ _       
                           |    /|\           /|\          |   /|\  
                           |     |penetration  |R_fuel     |    |
                           |     |             |           |    |
                   ________|_____|_____________|___________|____|__________
                  /|\      |                   |           |    |relocation*R_fuel_cold
                   |       |                   |           |    |
                   |Ri_clad|                   |        _ _|_ _\|/_ _ _ _ _ 
                   |       |                   |           |  
                */

                // Compute relocation recovery as a fraction of relocation
                scalar newRelRecover = relocation*(recFrac_);

                // Perturb epsilon Relocation for sensitivity analysis
                newRelRecover = newRelRecover * F_relocationRecovery.value() 
                            + delta_relocationRecovery.value();
                
                // Relocation is recovered as much as possible to prevent penetration.
                // If there is no penetration, there is no recovery
                newRelRecover = min
                (
                    newRelRecover, 
                    2/DiamCold_*penetration + relRecoverOld[cellI]
                );

                // Relax relocation recovery
                relRecover[cellI] = relaxRecovery_* newRelRecover 
                                + (1 - relaxRecovery_)*relRecover[cellI];

                rel[cellI] = 
                max
                (
                    relocation*(1 - recFrac_), 
                    (relocation - relRecover[cellI])
                );

            }

        }
    }
}


Foam::scalar Foam::relocationFrapcon::setRelocation
(
        double Q,
        double Gap,
        double D,
        double Buslice
)
{    
    double Bu = Buslice/1000/0.881;
    
    // Power in kW/m
    scalar Qprime = Q*constant::mathematical::pi*pow(D/2.0,2.0)/1000; 
    
    // Initialize gap change    
    scalar deltaGap(0.0);
    
    if(modifiedRelocationModel_)
    {
        if(Bu < 0.0937) // TODO : correction at low burnup or not?
        {
                deltaGap = 0.055*100;
        }
        else
        {
            // Inizialize reloc variable
            scalar reloc = 0;
            
            if(Qprime < 20)
            {
                    reloc = 0.345;
            }
            else if(Qprime <= 40)
            {
                    reloc = 0.345 + (Qprime - 20.0)/200.0;            
            }
            else
            {
                    reloc = 0.445;              
            }
        
            deltaGap = 
            (0.055 + min(reloc, reloc*(0.5795 + 0.2447*Foam::log(Bu))))*100;
        }
    }
    else
    {    
        scalar PFACTOR = (Qprime - 20)*5/20;
        
        scalar FBu = 0;
        if(Bu < 5)
        {
                FBu = Bu/5;
        }
        else
        {
                FBu = 1;
        }        

        if(Qprime < 20 + 1e-6 )
        {
                deltaGap = (30 + 10*FBu);            
        }
        else if(Qprime < 40 + 1e-6 )
        {
                deltaGap = (28 + PFACTOR + (12 + PFACTOR)*FBu);            
        }
        else
        {
                deltaGap = (32 + 18*FBu);            
        }
    }   
    
    const scalar nominalValue = deltaGap*(Gap/D)/100;

    return nominalValue;
     
}


// ************************************************************************* //
