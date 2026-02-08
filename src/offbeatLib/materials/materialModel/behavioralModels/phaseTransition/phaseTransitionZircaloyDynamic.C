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

#include "phaseTransitionZircaloyDynamic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseTransitionZircaloyDynamic, 0);
    addToRunTimeSelectionTable(phaseTransitionModel, phaseTransitionZircaloyDynamic, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransitionZircaloyDynamic::phaseTransitionZircaloyDynamic
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    phaseTransitionModel(mesh, dict, baseDict),
    w_(0.0)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransitionZircaloyDynamic::~phaseTransitionZircaloyDynamic()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseTransitionZircaloyDynamic::correct
(
    const scalarField& T, 
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_betaFraction") ?
    this->F_betaFraction() : 
    phaseTransitionModel::F_betaFraction();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_betaFraction") ?
    this->delta_betaFraction() : 
    phaseTransitionModel::delta_betaFraction();

    // Reference to oldTime and current betaFraction internalField
    const scalarField& betaFracOldI = betaFraction_.oldTime().internalField();
    scalarField&       betaFracI    = betaFraction_.ref();

    // Timestep value
    const scalar dt = mesh_.time().deltaT().value(); 

    // Old temperature field 
    const scalarField& TOld = 
        mesh_.lookupObject<volScalarField>("T").oldTime().internalField();

    // Compute central and span temperatures for the phase transition curve
    const scalar Tcent = 1159 - 0.096 * w_ ; 
    const scalar Tspan =   44 + 0.026 * w_ ;

    scalar relax = mesh_.relaxField(betaFraction_.name()) ? 
    mesh_.fieldRelaxationFactor(betaFraction_.name()) : 1.0;

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const scalar   Ti = T[cellI];

        // Compute rate of change of T
        const scalar Q = (Ti - TOld[cellI])/dt;

        // Compute onset temperature alpha->beta
        scalar Tab(-GREAT);
        if ( Q >= 0.0 and Q < 0.1 )
        {
            Tab = 1083 - 0.152 * w_;
        }
        else if ( Q >= 0.1 and Q <= 100 )
        {
            Tab = (1113 - 0.156*w_) * pow(Q,0.0118);
        }
        else if ( Q > 100 )
        {
            // Q is limited to 100 (out of correlation bounds)
            Tab = (1113 - 0.156*w_) * pow(100,0.0118);            
        }

        // Compute onset temperature beta->alpha
        scalar Tba(GREAT);
        if ( Q <= 0 and Q > -0.1 )
        {
            Tba = 1300;
        }
        else if ( Q <= -0.1 and Q >= -100 )
        {
            Tba = 1302.8 - 8.333 * pow(mag(Q),0.477);
        }
        else if ( Q < -100 )
        {
            // Q is limited to -100 (out of correlation bounds)
            Tba = 1302.8 - 8.333 * pow(mag(-100),0.477);
        }
        
        if( mag(Q) > 100 )
        {
            // Warning for validity bounds correlation
            WarningIn("void Foam::phaseTransitionZircaloyDynamic::correct") << nl  
            << "    Rate of change of temperature overcomes the validity bounds" 
            << " for this correlation."             << nl
            << " Max. acceptable rate : 100 K/s"    << nl
            << " Current rate : " << Q << " K/s"    << nl << endl;
        }

        if ( Ti > Tba )
        {
            betaFracI[cellI] = relax*1 + (1 - relax)*betaFracI[cellI];
        }
        else if ( Ti < Tab )
        {
            betaFracI[cellI] = relax*0 + (1 - relax)*betaFracI[cellI];
        }
        else
        {
            // Compute phase equilibrium fraction ys, Tcent, Tspan
            const scalar    ys = 0.5 * (1 + tanh((Ti-Tcent)/Tspan)); 

            // k0, km and k parameter
            const scalar km = ( Q > 0 ) ? 0.0 : 0.2;
            const scalar k0 = 60457 + 18129 * mag(Q);
            const scalar k  = (k0 * exp(-16650/Ti) + km);

            // Solve phase transformation ODE
            // Numerical integration with second order Adams-Moulton method (AM2)
            const scalar yOld   = betaFracOldI[cellI]; 
            scalar yNew   = pow(1+dt*k/2, -1) * (dt*k*ys + yOld*(1-dt*k/2));

            // Perturb epsilon densification for sensitivity analysis
            yNew = yNew * F.value() + delta.value();

            betaFracI[cellI] = relax*yNew + (1 - relax)*betaFracI[cellI];
        }

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            {
                if (patchID > -1 and betaFraction_.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                    scalarField& betaFracP = 
                    betaFraction_.boundaryFieldRef()[patchID];

                    const scalarField& betaFracOldP = 
                    betaFraction_.oldTime().boundaryField()[patchID];

                    const scalarField& Tp = 
                    mesh_.lookupObject<volScalarField>("T").boundaryField()[patchID];

                    const scalarField& TOldP = 
                    mesh_.lookupObject<volScalarField>("T").oldTime().boundaryField()[patchID];

                    const scalar   Ti = Tp[faceID];

                    //- Compute rate of change of T
                    const scalar Q = (Ti - TOldP[faceID])/mesh_.time().deltaT().value();

                    // Compute onset temperature alpha->beta
                    scalar Tab(-GREAT);
                    if ( Q >= 0.0 and Q < 0.1 )
                    {
                        Tab = 1083 - 0.152 * w_;
                    }
                    else if ( Q >= 0.1 and Q <= 100 )
                    {
                        Tab = (1113 - 0.156*w_) * pow(Q,0.0118);
                    }
                    else
                    {
                        // Q is limited to 100 (out of correlation bounds)
                        Tab = (1113 - 0.156*w_) * pow(100,0.0118);            
                    }

                    // Compute onset temperature beta->alpha
                    scalar Tba(GREAT);
                    if ( Q > -0.1 and Q <= 0 )
                    {
                        Tba = 1300;
                    }
                    else if ( Q >= -100 and Q <= -0.1 )
                    {
                        Tba = 1302.8 - 8.333 * pow(mag(Q),0.477);
                    }
                    else
                    {
                        // Q is limited to -100 (out of correlation bounds)
                        Tba = 1302.8 - 8.333 * pow(mag(-100),0.477);
                    }

                    if ( Ti > Tba )
                    {
                        betaFracP[faceID] = relax*1 + (1 - relax)*betaFracP[faceID];
                    }
                    else if ( Ti < Tab )
                    {
                        betaFracP[faceID] = relax*0 + (1 - relax)*betaFracP[faceID];
                    }
                    else
                    {
                        // Compute phase equilibrium fraction ys 
                        const scalar    ys = 0.5 * (1 + tanh((Ti-Tcent)/Tspan)); 

                        // k0, km and k parameter
                        const scalar km = ( Q > 0 ) ? 0.0 : 0.2;
                        const scalar k0 = 60457 + 18129 * mag(Q);
                        const scalar k  = (k0 * exp(-16650/Ti) + km);

                        // Solve phase transformation ODE
                        // Numerical integration with second order Adams-Moulton method (AM2)
                        const scalar yOld   = betaFracOldP[faceID]; 
                        scalar yNew   = pow(1+dt*k/2, -1) * (dt*k*ys + yOld*(1-dt*k/2));
                        
                        // Perturb epsilon densification for sensitivity analysis
                        yNew = yNew * F.value() + delta.value();
                        
                        betaFracP[faceID] = relax*yNew + (1 - relax)*betaFracP[faceID];
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
