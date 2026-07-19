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

#include "PuRedistributionSolver.H"
#include "zeroCurrentActinidesRedistributionFvPatchScalarField.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PuRedistributionSolver, 0);
    addToRunTimeSelectionTable
    (
        PuRedistributionSolver, 
        PuRedistributionSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PuRedistributionSolver::initializePuConcentration()
{
    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const fuelMaterial& fuelMat = 
                refCast<const fuelMaterial>(mat_.materialsList()[i]);
            
            // Lookup initial Pu weight fraction
            const dictionary& PuDict = fuelMat.isotopesDict().subDict("Pu");
            const scalar Pu0(readScalar(PuDict.lookup("ratioOverMetal")));
            
            // Assign to each fuel cell
            const labelList& addr(matAddrList[i]);
            forAll(addr, j)
            {
                const label cellI = addr[j];
                PuFormFactor_[cellI] = Pu0;
            }
        }
    }

    PuFormFactor_.correctBoundaryConditions();
}

void Foam::PuRedistributionSolver::updateCoefficients()
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
                const scalar OMth    = 2.0-pow((16.0*exp(- min(3.3/KbT, 75))),0.25);
                const scalar Nld194  = 0.0009*exp(- min(3.1/KbT, 75));
                const scalar NldOMth = (4.0/pow((2.0-OMth),2.0))*exp(-min(6.4/KbT,75));
                const scalar Nld200  = 2.0*exp(-min(4.2/KbT,75));

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

                // if( Ti[cellI] > 2173 )
                // {
                    // diffPu_[cellI] =
                    //     fOM * 3.4e-5 * 3 * exp(-min(55891/Ti[cellI], 50)); 
                // }
                // else
                // {
                //     diffPu_[cellI] = 
                //         fOM * 1.5e-10 * exp(-min(2.22e5/(8.31446*Ti[cellI]), 50));
                // }

                fOM = min(fOM, 1.0);
                fOM = max(fOM, 1e-5);

                diffPu_[cellI] =
                    fOM * 3.4e-5 * exp(-min(55891/Ti[cellI], 50)); 
                


                // --- FACE VALUES --- //

                // Loop over all faces of current cell
                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    
                    if (patchID > -1 and PuFormFactor_.boundaryField()[patchID].size())
                    {
                        // Take references to patch fields
                        scalarField& diffPuP(diffPu_.boundaryFieldRef()[patchID]);
                        const scalarField& Tp(T_.boundaryField()[patchID]);

                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                        // Parameters to compute correction factor fOM
                        const scalar KbT     = 8.617e-5*Tp[faceID];
                        const scalar OMth    = 2.0-pow((16.0*exp(- min(3.3/KbT, 75))),0.25);
                        const scalar Nld194  = 0.0009*exp(-min(3.1/KbT,75));
                        const scalar NldOMth = (4.0/pow((2.0-OMth),2.0))*exp(-min(6.4/KbT,75));
                        const scalar Nld200  = 2.0*exp(-min(4.2/KbT,75));

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

                        fOM = min(fOM, 1.0);
                        fOM = max(fOM, 1e-5);

                        // if( Tp[faceID] > 2173 )
                        // {
                            // diffPuP[faceID] =
                                // fOM * 3.4e-5 * 3 * exp(-min(55891/Tp[faceID],50)); 
                        // }
                        // else
                        // {
                        //     diffPuP[faceID] = 
                        //         fOM * 1.5e-10 * exp(-min(2.22e5/(8.31446*Tp[faceID]),50));
                        // }
                        
                        diffPuP[faceID] =
                            fOM * 3.4e-5 * exp(-min(55891/Tp[faceID],50)); 
                    }
                }
            }
        }
    }

    diffPu_.correctBoundaryConditions();
}

void Foam::PuRedistributionSolver::redistributeIsotopes()
{
    volScalarField& N_Pu(mesh_.lookupObjectRef<volScalarField>("N_Pu"));
    volScalarField& N_Pu238(mesh_.lookupObjectRef<volScalarField>("N_Pu238"));
    volScalarField& N_Pu239(mesh_.lookupObjectRef<volScalarField>("N_Pu239"));
    volScalarField& N_Pu240(mesh_.lookupObjectRef<volScalarField>("N_Pu240"));
    volScalarField& N_Pu241(mesh_.lookupObjectRef<volScalarField>("N_Pu241"));
    volScalarField& N_Pu242(mesh_.lookupObjectRef<volScalarField>("N_Pu242"));

    const labelListList& matAddrList(mat_.matAddrList());
    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                const scalar ratio = 
                    PuFormFactor_[cellI] / PuFormFactor_.oldTime()[cellI];

                // Pu Total
                N_Pu[cellI] = N_Pu.oldTime()[cellI] * ratio;

                N_Pu238[cellI] = N_Pu238.oldTime()[cellI] * ratio;
                N_Pu239[cellI] = N_Pu239.oldTime()[cellI] * ratio;
                N_Pu240[cellI] = N_Pu240.oldTime()[cellI] * ratio;
                N_Pu241[cellI] = N_Pu241.oldTime()[cellI] * ratio;
                N_Pu242[cellI] = N_Pu242.oldTime()[cellI] * ratio;
            }
        }
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PuRedistributionSolver::PuRedistributionSolver
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
    PuFormFactor_
    (
        IOobject
        (
            "PuFormFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        "zeroGradient"
    ),
    diffPu_
    (
        IOobject
        (
            "diffPu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
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
    // Initialize Pu concentration (if not already read from field, e.g. restart
    // simulation)
    if (max(PuFormFactor_.internalField()).value() == 0)
    {
        initializePuConcentration();
    }

    // Correct coefficients
    updateCoefficients();

    // Read options PuRedistribution from dedicated subdict 
    dictionary PuDict
    (
        elementTransportDict.
        subOrEmptyDict("PuRedistributionOptions")
    );        
    
    poreMigrationOn_ = 
    PuDict.lookupOrDefault<bool>
    (
        "poreMigration", false
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::PuRedistributionSolver::correct()
{
    // Update thermal diffusion coefficient
    updateCoefficients();

    const int nCorr = transportSolver::nCorrectors(PuFormFactor_.name());

    int nInnerIter = 0;

    bool convergedInner(false);

    do
    {
        // Store prev iter values for rel residual
        PuFormFactor_.storePrevIter();

        // Flux for Soret term
        surfaceScalarField soretTerm = 
        // fvc::interpolate( diffPu_*Q_/(R_*pow(T_,2.0)) * gradT_ ) &  mesh_.Sf();
        fvc::interpolate( diffPu_*Q_/(R_*pow(T_,2.0)) * gradT_ * (1-PuFormFactor_)) &  mesh_.Sf();  

        // Equation for Pu redistribution
        fvScalarMatrix PuEqn
        (
              fvm::ddt(PuFormFactor_)
            - fvm::laplacian(diffPu_, PuFormFactor_)
            - fvm::div( soretTerm, PuFormFactor_ , "div(soretTerm,PuFormFactor)")
            // + fvc::div( soretTerm, pow(PuFormFactor_, 2.0) )
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
                porosity_ * diffPu_ / poreDiameter_ * A_ * poreThickness_
                * Foam::exp(-diffPu_ / (poreThickness_*max(mag(pV),verySmall)))
                * gradT_ 
            ) & mesh_.Sf();
            
            PuEqn += fvm::div(pmTerm, PuFormFactor_, "div(pmTerm,PuFormFactor)");
        }    
        
        // Relax equation    
        PuEqn.relax();
        
        // Compute residual
        residual_ = PuEqn.solve().max().initialResidual();

        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
        }
        
        // Relax fields
        PuFormFactor_.relax();

        // Calculate a different residual based on the relative change
        scalar denom = 
        gMax(mag(PuFormFactor_.primitiveField() - PuFormFactor_.oldTime().primitiveField()));

        scalar num = 
        gMax(mag(PuFormFactor_.primitiveField() - PuFormFactor_.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(PuFormFactor_.primitiveField())), SMALL);
        }
        
        relResidual_ = num/denom;

        if(denom < 1e-5)
        {
            relResidual_ = VSMALL;
        }

        Info << "relResidualPu " << relResidual_ << endl;

        convergedInner = converged(PuFormFactor_.name());

    } while 
    (
        not(convergedInner)  
        && ++nInnerIter < nCorr 
    );

    // Redistribute all the isotopes + total Pu
    redistributeIsotopes();

    Info << endl;
}

// ************************************************************************* //
