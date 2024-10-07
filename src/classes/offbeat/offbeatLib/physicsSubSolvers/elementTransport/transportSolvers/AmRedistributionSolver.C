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

#include "AmRedistributionSolver.H"
#include "zeroCurrentActinidesRedistributionFvPatchScalarField.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(AmRedistributionSolver, 0);
    addToRunTimeSelectionTable
    (
        AmRedistributionSolver, 
        AmRedistributionSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::AmRedistributionSolver::initializeAmConcentration()
{
    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const fuelMaterial& fuelMat = 
                refCast<const fuelMaterial>(mat_.materialsList()[i]);
            
            const scalar Am0 = fuelMat.ratioAmMetal();

            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                AmFormFactor_[cellI] = Am0;
                
            }
        }
    }

    AmFormFactor_.correctBoundaryConditions();
}

void Foam::AmRedistributionSolver::updateCoefficients()
{
    const scalarField& Ti(T_.internalField());

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const fuelMaterial& fuelMat = 
                refCast<const fuelMaterial>(mat_.materialsList()[i]);
            
            const scalar OM = fuelMat.oxygenMetalRatio();

            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                // --- Diffusion Coefficient --- //

                // Parameters to compute correction factor fOM
                const scalar KbT     = 8.617e-5*Ti[cellI];
                const scalar OMth    = 2.0-pow((16.0*exp(-3.3/KbT)),0.25);
                const scalar Nld194  = 0.0009*exp(-3.1/KbT);
                const scalar NldOMth = (4.0/pow((2.0-OMth),2.0))*exp(-6.4/KbT);
                const scalar Nld200  = 2.0*exp(-4.2/KbT);

                // Compute correcting factor fOM
                const scalar fOM194 = Nld194 / Nld200; 
                const scalar fOMth  = NldOMth / Nld200;

                scalar fOM(-GREAT);
                
                if( OM <= 1.94 )
                {
                    fOM = fOM194;
                }
                else if( OM < OMth ) 
                {
                    // Log interpolation
                    fOM = exp
                    (
                        log(fOM194) 
                        + (log(OM)-log(1.94))/(log(OMth)-log(1.94))
                        * (log(fOMth)-log(fOM194))
                    );
                }
                else if( OM == OMth )
                {
                    fOM = fOMth;
                }
                else if( OM < 2.00 )
                {
                    // Log interpolation
                    fOM = exp
                    (
                        log(fOMth) 
                        + (log(OM)-log(OMth))/(log(2.0)-log(OMth))
                        * (log(1.0)-log(fOMth))
                    );
                }
                else if( OM >= 2.00 )
                {
                    fOM = 1.00;
                }

                diffAm_[cellI] = fOM * 2.0e-6 * exp(-46324/Ti[cellI]);

                // --- FACE VALUES --- //

                // Loop over all faces of current cell
                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    
                    if (patchID > -1 and AmFormFactor_.boundaryField()[patchID].size())
                    {
                        // Take references to patch fields
                        scalarField& diffAmP(diffAm_.boundaryFieldRef()[patchID]);
                        const scalarField& Tp(T_.boundaryField()[patchID]);

                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                        // Parameters to compute correction factor fOM
                        const scalar KbT     = 8.617e-5*Tp[faceID];
                        const scalar OMth    = 2.0-pow((16.0*exp(-3.3/KbT)),0.25);
                        const scalar Nld194  = 0.0009*exp(-3.1/KbT);
                        const scalar NldOMth = (4.0/pow((2.0-OMth),2.0))*exp(-6.4/KbT);
                        const scalar Nld200  = 2.0*exp(-4.2/KbT);

                        // Compute correcting factor fOM
                        const scalar fOM194 = Nld194 / Nld200; 
                        const scalar fOMth  = NldOMth / Nld200;

                        scalar fOM(-GREAT);
                        
                        if( OM <= 1.94 )
                        {
                            fOM = fOM194;
                        }
                        else if( OM < OMth ) 
                        {
                            // Log interpolation
                            fOM = exp
                            (
                                log(fOM194) 
                                + (log(OM)-log(1.94))/(log(OMth)-log(1.94))
                                * (log(fOMth)-log(fOM194))
                            );
                        }
                        else if( OM == OMth )
                        {
                            fOM = fOMth;
                        }
                        else if( OM < 2.00 )
                        {
                            // Log interpolation
                            fOM = exp
                            (
                                log(fOMth) 
                                + (log(OM)-log(OMth))/(log(2.0)-log(OMth))
                                * (log(1.0)-log(fOMth))
                            );
                        }
                        else if( OM >= 2.00 )
                        {
                            fOM = 1.00;
                        }

                        diffAmP[faceID] =
                            fOM * 2.0e-6 * exp(-46324/Tp[faceID]);

                    }
                }
            }
        }
    }

    diffAm_.correctBoundaryConditions();
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AmRedistributionSolver::AmRedistributionSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    transportSolver(mesh, mat, elementTransportDict, solverName),
    T_(mesh_.lookupObject<volScalarField>("T")),
    gradT_(mesh_.lookupObject<volVectorField>("gradT")),
    AmFormFactor_
    (
        IOobject
        (
            "AmFormFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        "zeroGradient"
        // "zeroCurrentActinidesRedistribution"
    ),
    diffAm_
    (
        IOobject
        (
            "diffAm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, 0.0),
        "calculated"
    ),
    Q_("Q", dimEnergy/dimMoles, -146.5e3),
    R_("R", dimEnergy/dimMoles/dimTemperature, 8.31446),
    porosity_(mesh_.lookupObject<volScalarField>("porosity")),
    poreDiameter_("poreDiameter", dimLength, 0.008e-3),
    poreThickness_("poreThickness", dimLength, 80e-3),
    A_("A", dimless/dimTemperature, 0.35),
    poreMigrationOn_()
{

    // Initialize Am concentration
    initializeAmConcentration();

    // Correct coefficients
    updateCoefficients();

    // Read options AmRedistribution from dedicated subdict 
    dictionary AmDict
    (
        elementTransportDict.
        subOrEmptyDict("AmRedistributionOptions")
    );        
    
    poreMigrationOn_ = 
    AmDict.lookupOrDefault<bool>
    (
        "poreMigration", true
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::AmRedistributionSolver::correct()
{
    // Update thermal diffusion coefficient
    updateCoefficients();

    const dictionary& stressControl =  mesh_.solutionDict().subDict("stressAnalysis");
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    int nInnerIter = 0;

    bool convergedInner(false);

    do
    {
        // Store prev iter values for rel residual
        AmFormFactor_.storePrevIter();

        // Flux for Soret term
        surfaceScalarField soretTerm = 
        fvc::interpolate( diffAm_*Q_/(R_*pow(T_,2.0)) * gradT_ ) &  mesh_.Sf();

        // Equation forAm redistribution
        fvScalarMatrix AmEqn
        (
              fvm::ddt(AmFormFactor_)
            - fvm::laplacian(diffAm_, AmFormFactor_)
            - fvm::div( soretTerm, AmFormFactor_ , "div(soretTerm,AmFormFactor)")
            + fvc::div( soretTerm, pow(AmFormFactor_, 2.0) )
        );

        if (poreMigrationOn_)
        {
            // Useful dimensioned scalar
            const dimensionedScalar verySmall("verySmall", dimLength/dimTime, VSMALL);

            // Take reference to pore velocity
            const volVectorField& pV(mesh_.lookupObject<volVectorField>("poreVelocity"));

            // Flux for pore migration term
            surfaceScalarField pmTerm = 
            fvc::interpolate
            (
                porosity_ * diffAm_ / poreDiameter_ * A_ * poreThickness_
                * Foam::exp(-diffAm_ / (poreThickness_*max(mag(pV),verySmall)))
                * gradT_ 
            ) & mesh_.Sf();

            AmEqn += fvm::div(pmTerm, AmFormFactor_, "div(pmTerm,AmFormFactor)");
        }

        // Relax equation    
        AmEqn.relax();
        
        // Compute residual
        residual_ = AmEqn.solve().max().initialResidual();

        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
        }
        
        // Relax fields
        AmFormFactor_.relax();

        // Calculate a different residual based on the relative change
        scalar denom = 
        gMax(mag(AmFormFactor_.primitiveField() - AmFormFactor_.oldTime().primitiveField()));

        scalar num = 
        gMax(mag(AmFormFactor_.primitiveField() - AmFormFactor_.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(AmFormFactor_.primitiveField())), SMALL);
        }
        
        relResidual_ = num/denom;

        if(denom < 1e-5)
        {
            relResidual_ = VSMALL;
        }

        Info << "relResidualAm " << relResidual_ << endl;

        convergedInner = converged(AmFormFactor_.name());

    } while 
    (
        not(convergedInner)  
        && ++nInnerIter < nCorr 
    );

    Info << endl;
}

// ************************************************************************* //
