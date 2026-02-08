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

#include "swellingFBRMOX.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingFBRMOX, 0);
    addToRunTimeSelectionTable(swellingModel, swellingFBRMOX, dictionary);

    const char* swellingFBRMOX::group_ = 
        ("behavioralModels::swelling::FBRMOX");

    defineParameter(swellingFBRMOX, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingFBRMOX, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::scalarField Foam::swellingFBRMOX::computeGapWidth1D()
{
    scalarField gw1D(mapper_->nSlices(), 0.0);
    
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

            List<cell> sliceCells(mesh_.cells(), sliceAddr);
            scalarField sliceGapWidth(sliceAddr.size(), 0.0);
            scalarField sliceMagSf(sliceAddr.size(), 0.0);

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

                        sliceGapWidth[j] += gapWidthP[faceID]*faceSf;
                        sliceMagSf[j] += faceSf;
                    }  
                }
            }

            gw1D[i] = gSum(sliceGapWidth*sliceMagSf)/gSum(sliceMagSf);
        }
    }
    return gw1D;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingFBRMOX::swellingFBRMOX
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    swellingModel(mesh, dict, baseDict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    gapWidthName_(dict.lookupOrDefault<word>("gapWidthName", "gapWidth" )),
    gapWidth_(nullptr),
    outerPatchID_(-1),
    par1(0.020),
    par2(0.012),
    par3(0.0065),
    mapper_(nullptr),
    intraGasSwell_(nullptr),
    interGasSwell_(nullptr)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 0.020);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", 0.012);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", 0.0065);
    }

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

Foam::swellingFBRMOX::~swellingFBRMOX()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingFBRMOX::correct
(
    const scalarField& T,
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F =
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() :
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta =
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() :
    swellingModel::delta_epsilonSwelling();

    const fvMesh& mesh_ = epsilonSwelling_.mesh();

    if( Bu_ == nullptr 
        or mapper_ == nullptr 
        or gapWidth_ == nullptr
        or intraGasSwell_ == nullptr
        or interGasSwell_ == nullptr
    )
    {
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
        gapWidth_ = &mesh_.lookupObject<volScalarField>(gapWidthName_);
        mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
        if(mesh_.foundObject<volSymmTensorField>("intragranularGasSwelling"))
        {
            intraGasSwell_ = &mesh_.lookupObject<volSymmTensorField>("intragranularGasSwelling");
            interGasSwell_ = &mesh_.lookupObject<volSymmTensorField>("intergranularGasSwelling");
        }
    }

    // Get n of slices
    scalar nSlices(mapper_->nSlices());
    if( nSlices <= 0 )
    {
        FatalErrorInFunction()
        << "Number of slices cannot be less than 1." << nl
        << "Check that a sliceMapper different than \"none\" is used "
        << "or check the slice definition."
        << abort(FatalError);
    }

    // Take const references Bu fields
    const volScalarField& Bu = mesh_.lookupObject<volScalarField>(burnupName_);

    // Calculate slice-average values
    const scalarField& Bu_1D(mapper_->sliceAverage(Bu));
    const scalarField& BuOld_1D(mapper_->sliceAverage(Bu.oldTime()));

    // Compute gapWidth slice-wise field (each slice the fuel outer patch value)
    scalarField gapWidth_1D = computeGapWidth1D();

    // Ref to swelling fields
    const symmTensorField& swellingIOld = epsilonSwelling_.oldTime().internalField();
    symmTensorField& swellingI = epsilonSwelling_.ref();

    // Slice IDs field
    const scalarField& sliceIDs(mapper_->sliceID().internalField());

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const label sliceID = sliceIDs[cellI];

        bool gapOpen = ( gapWidth_1D[sliceID] <= 0 )
        ? false
        : true;

        // Convert Bu in MWd/tUO2 to FIMA %
        const scalar BuI = Bu_1D[sliceID] * 1e-3 / 9.5;
        const scalar BuOldI = BuOld_1D[sliceID] * 1e-3 / 9.5;

        // Compute Bu increment
        const scalar deltaBu = BuI-BuOldI;

        // Compute solid swelling
        const scalar solidSwellIncrement = par3*deltaBu;
        
        // Init total swelling increment
        scalar totSwellIncrement = 0;

        if ( BuI <= 1 and gapOpen )
        {
            totSwellIncrement = par1*deltaBu;

        }
        else if ( BuI > 1 and gapOpen )
        {
            totSwellIncrement = par2*deltaBu;

        }
        else // gap closed
        {
            totSwellIncrement = solidSwellIncrement;
        }

        // Perturb
        totSwellIncrement = totSwellIncrement * F.value() + delta.value();

        // if(mesh_.foundObject<volSymmTensorField>("intragranularGasSwelling"))
        // {
        //     // Assign total swelling
        //     swellingI[cellI] = 
        //     swellingIOld[cellI] + 
        //     (
        //         solidSwellIncrement/3.0 * I +
        //         (intraGasSwell_->internalField()[cellI]-intraGasSwell_->oldTime().internalField()[cellI])/3 +
        //         (interGasSwell_->internalField()[cellI]-interGasSwell_->oldTime().internalField()[cellI])/3
        //     );
        // }
        // else
        // {
            // Assign total swelling
            swellingI[cellI] = 
                swellingIOld[cellI] + (totSwellIncrement)/3.0 * I;
        // }

    }
}


// ************************************************************************* //
