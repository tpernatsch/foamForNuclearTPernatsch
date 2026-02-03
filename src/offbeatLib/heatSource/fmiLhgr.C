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

#if defined __has_include
#  if __has_include(<commDataLayer.H>) 
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "commDataLayer.H"
#include "fmiLhgr.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "IndirectList.H"
    
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fmiLhgr, 0);
    
    addToRunTimeSelectionTable
    (
        heatSource, 
        fmiLhgr, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fmiLhgr::fmiLhgr
(
    const fvMesh& mesh,
    const materials& materials,
    const dictionary& heatSourceOptDict,
    const bool timeDependent
)
:
    constantLhgr(mesh, materials, heatSourceOptDict, true),
    fmiState_
    (
        IOobject
        (
            "lhgr",
            mesh_.time().timeName(),
            "uniform/fmiState",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    lhgrAvg_(fmiState_.lookupOrDefault("lhgrAvg", 0.0)),
    Qavg_(fmiState_.lookupOrDefault("Qavg", 0.0)),
    alpha_(heatSourceOptDict.lookupOrDefault("relationFactor", 1.0)),
    isUseLhgr_(false),
    isUseQavg_(false)
{
    // Preparing data to update temperature
    word lhgrKeyFromFMU("lhgrNameFromFMU");
    word QavgKeyFromFMU("QavgNameFromFMU");
    if (heatSourceOptDict.found(lhgrKeyFromFMU))
    {
        isUseLhgr_ = true;

        lhgrNameFromFMU_ = heatSourceOptDict.get<word>(lhgrKeyFromFMU);
        
        // Communicating with the FMU
        const Time& runTime = mesh_.time();
        commDataLayer& data = commDataLayer::New(runTime); 

        // Store in data layer and set its initial value to the T 
        // in the dictionary 
        data.storeObj(
            lhgrAvg_,
            lhgrNameFromFMU_,
            commDataLayer::causality::in
        );
    }
    else if (heatSourceOptDict.found(QavgKeyFromFMU))
    {
        isUseQavg_ = true;

        QavgNameFromFMU_ = heatSourceOptDict.get<word>(QavgKeyFromFMU);
        
        // Communicating with the FMU
        const Time& runTime = mesh_.time();
        commDataLayer& data = commDataLayer::New(runTime); 

        // Store in data layer and set its initial value to the T 
        // in the dictionary 
        data.storeObj(
            Qavg_,
            QavgNameFromFMU_,
            commDataLayer::causality::in
        );
    }
    else
    {
        FatalIOErrorInFunction(heatSourceOptDict)
            << "Cannot find 'lhgrNameFromFMU' or 'QavgNameFromFMU' key in "
            << "heatSourceOptions"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fmiLhgr::~fmiLhgr()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fmiLhgr::correct()
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();

    if(mesh_.topoChanging())
    {
        calcAddressing(heatSourceOptDict_);
        calcReferenceDimensions(heatSourceOptDict_);
    }

    if(!mesh_.foundObject<volScalarField>("QPrevIter"))
    {
        Q_.storePrevIter();
    }

    Q_.prevIter().storePrevIter();        
    Q_.storePrevIter();

    const scalarField& Qprev(Q_.prevIter());
    const scalarField& QprevPrev(Q_.prevIter().prevIter());

    // Correct heat source only once per time step or until power changes
    if
    (
        gMax(mag(QprevPrev - Qprev)) > 1e-3 
        or
        mesh_.time().value() > currentTime_
        or
        mesh_.foundObject<volScalarField>("porosity")
    )
    {  
        currentTime_ = mesh_.time().value();

        scalarField Vi(mesh_.V(), addr_);
        scalar totalVol(gSum(Vi)/angularFraction);

        // Communicating with the FMU
        const Time& runTime = mesh_.time();
        commDataLayer& data = commDataLayer::New(runTime);

        // Get new_step flag, if true apply direcly the LHGR value from the FMI, 
        // else apply the relaxation factor
        label isNewStep = data.getObj<label>
        (
            "new_step", commDataLayer::causality::in
        );
        
        if (isUseQavg_)
        {
            scalar QavgFMI = data.getObj<scalar>
            (
                QavgNameFromFMU_, 
                commDataLayer::causality::in
            );

            Qavg_ = isNewStep
                ? QavgFMI
                : alpha_*QavgFMI + (1.0-alpha_)*Qavg_;

            lhgrAvg_ = Qavg_*totalVol/(zMax_ - zMin_);
        }
        else if (isUseLhgr_)
        {
            scalar lhgrFMI = data.getObj<scalar>
            (
                lhgrNameFromFMU_, 
                commDataLayer::causality::in
            );

            lhgrAvg_ = isNewStep
                ? lhgrFMI
                : alpha_*lhgrFMI + (1.0-alpha_)*lhgrAvg_;

            // Interpolate LHGR and calculate average power density Qavg
            Qavg_ = lhgrAvg_*(zMax_ - zMin_)/totalVol;
        }

        // Correct the radial and axial profiles
        radialModel_->correct();
        axialModel_->correct();

        // Build the resulting heat source
        scalarField& Qi = Q_.ref();
        scalarField& lhgri = lhgr_.ref();
        const scalarField& f_r = radialModel_->profile();
        const scalarField& g_z = axialModel_->profile();

        forAll(addr_, addrI)
        {
            const label cellI = addr_[addrI];
            
            Qi[cellI] = Qavg_*f_r[addrI]*g_z[addrI];
            lhgri[cellI] = lhgrAvg_*g_z[addrI];
        }

        // Reshape power profile due to porosity redistribution
        reshapePowerPorosity();
        
        scalarField Qaddr(Q_.ref(), addr_);
        scalar powerTot   = lhgrAvg_*(zMax_ - zMin_)*angularFraction;
        scalar powerTotOF = gSum( Vi*Qaddr );
        scalar error = 
        100*
        (
            mag(powerTotOF - powerTot)/max(powerTot, VSMALL)
        );

        if(error > 0.01)
        {            
            WarningIn("Foam::fmiLhgr::correct()") << nl  
                << "    Total power does not correspond to LHGR*height, "<< nl 
                << "    with an error of "  << error << "%" << nl
                << "    Radial or axial profiles might be not normalized." << nl
                << "    The field Q will be scaled to preserve total power."
                <<  nl << endl;

            radialModel_ -> checkNormalization();
            axialModel_ -> checkNormalization();

            scalar scalingFactor(powerTot/powerTotOF);

            forAll(addr_, addrI)
            {
                const label cellI = addr_[addrI];
                
                Qi[cellI] *= scalingFactor;
            }
        }
    }
    else
    {
        currentTime_ = mesh_.time().value();
    }

    Q_.correctBoundaryConditions();
    lhgr_.correctBoundaryConditions();

    //- Write FMI state
    if (mesh_.time().writeTime())
    {
        fmiState_.set("lhgrAvg", lhgrAvg_);
        fmiState_.set("Qavg", Qavg_);
    }
}


#endif

// ************************************************************************* //
