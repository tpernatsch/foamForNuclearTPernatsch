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

#include "HNGD.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"
#include "volFields.H"
#include "processorFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members  * * * * * * * * * * * * * //

defineTypeNameAndDebug(HNGD, 0);

addToRunTimeSelectionTable
(
    transportSolver,
    HNGD,
    dictionary
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Compute hydrides formation energy
scalar HNGD::hydrideFormationEnergy
(
    const scalar& T
)
{
    const scalar e = constant::electromagnetic::e.value();
    const scalar NA = constant::physicoChemical::NA.value();
    return (-0.5655 + 4e-4*T - 2e-7*pow(T,2) + 3e-10*pow(T,3)) * e * NA;
}


//- Convert concentration to atomic Fraction
scalar HNGD::convertToAtomFrac
(
    const scalar& c
)
{
     return c / (MH_ * (c/MH_ + (1e6 - c)/MZr_));
}


//- Compute volume fraction of alpha phase
scalar HNGD::fAlpha
(
    const scalar& T,
    const scalar& TSSd,
    const scalar& Cpp,
    const scalar& Ctot
)
{
    //- TSSd atomic fraction
    const scalar TSSd_at = convertToAtomFrac(TSSd);

    //- Hydride atomic fraction
    const scalar xhydride = convertToAtomFrac(Cpp);

    //- boundary fraction
    const scalar xdelta = hydrideBoundaryFraction(T);

    return 1 - min(xhydride /(xdelta - TSSd_at), 1); 
}


scalar HNGD::leverRule
(
    const scalar& T,
    const scalar& TSSd,
    const scalar& Ctot
)
{
    //- Total atomic fraction
    const scalar xtot = convertToAtomFrac(Ctot);

    //- TSSd atomic fraction
    const scalar xalpha = convertToAtomFrac(TSSd);

    //- alpha / alpha + delta boundary
    const scalar xdelta = hydrideBoundaryFraction(T);

    return ((xtot - xalpha) / (xdelta - xalpha));
}


//- Kinetic parameter for growth
scalar HNGD::KG
(
    const scalar& T,
    const scalar& TSSd,
    const scalar& Cpp,
    const scalar& Ctot
)
{
    const scalar e = constant::electromagnetic::e.value();
    const scalar NA = constant::physicoChemical::NA.value();
    const scalar R = constant::physicoChemical::R.value();
    
    const scalar Kmob0 = 5.53e5;
    const scalar Kth0 = 1.6e-5;
    const scalar Eg = 0.9 * e * NA;

    const scalar Kmob = Kmob0 * max(fAlpha(T, TSSd, Cpp, Ctot), SMALL)  
                        * leverRule(T, TSSd, Ctot) * exp(-Eg/ (R*T));

    const scalar Kth = Kth0 *  max(fAlpha(T, TSSd, Cpp, Ctot), SMALL) * leverRule(T, TSSd, Ctot) 
    * exp(-hydrideFormationEnergy(T) / (R*T));

    return pow((1/Kmob) + (1/Kth), -1);
}


//- Update material properties relative to the hydrogen transport 
void HNGD::updateCoefficients()
{
    // Constants
    const scalar e = constant::electromagnetic::e.value();
    const scalar NA = constant::physicoChemical::NA.value();
    const scalar R = constant::physicoChemical::R.value();
    const scalar Ed = 0.46 * e * NA;

    // Field references
    scalarField& KD = KD_.field();
    scalarField& KN = KN_.field();
    scalarField& KG = KG_.field();
    const scalarField& Cpp = Cpp_.oldTime().field();
    const scalarField& Css = Css_.oldTime().field();
    const scalarField& TSSd = TSSd_.field();
    const scalarField& T = T_.field();

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        if (TSSdModel_.set(i))
        {
            const labelList& addr(matAddrList[i]);

            //- Modif TSSd,TSSp + kinetic parameters calculations
            forAll(addr, i)
            {
                const label cellI = addr[i];

                // Update kinetics parameters
                KD[cellI] = 1.11e3 * exp(- Ed/ (R * T[cellI]));

                KN[cellI] = 2.75e-5 * fAlpha(T[cellI], TSSd[cellI], Cpp[cellI], (Cpp[cellI] + Css[cellI]))
                             * exp(-hydrideFormationEnergy(T[cellI])/(R * T[cellI]));

                KG[cellI] = HNGD::KG(T[cellI], TSSd[cellI], Cpp[cellI], (Cpp[cellI] + Css[cellI]));
            }
        }
    }

    KD_.correctBoundaryConditions();
    KN_.correctBoundaryConditions();
    KG_.correctBoundaryConditions();
}


void HNGD::updateRates()
{
    const scalarField& Cpp = Cpp_.oldTime().field();
    const scalarField& Css = Css_.oldTime().field();

    Kss_.field() = 0.0;
    S_.field() = 0.0;

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        if (TSSdModel_.set(i))
        {
            const labelList& addr(matAddrList[i]);
            forAll(addr, j)
            {
                const label cellI = addr[j];

                //- Nucleation Rate
                if (Css[cellI] >= TSSp_[cellI])
                {
                    Kss_[cellI] -= KN_[cellI];
                    S_[cellI] += KN_[cellI]*TSSp_[cellI];
                }
                //- Growth Rate
                else if ((Css[cellI] > TSSd_[cellI]) && Cpp[cellI] > SMALL)
                {
                    scalar x = Cpp[cellI] / max(Cpp[cellI] + Css[cellI] - TSSd_[cellI], SMALL);
                    
                    x = min(max(x, 0), 1-SMALL);

                    const scalar p = 2.5;
                    scalar K = KG_[cellI] * p * pow(-log(1-x), 1-1/p);

                    Kss_[cellI] -= K;
                    S_[cellI] += K*TSSd_[cellI];
                }
                //- Dissolution Rate
                //Comment/Uncomment depending on whether or not T condition is needed
                else if ((Css[cellI] <= TSSd_[cellI]) && (Cpp[cellI] > 0) /*&& (T_[cellI] - T0[cellI] > 0.01)*/)
                {
                    Kss_[cellI] -= KD_[cellI];
                    S_[cellI] += KD_[cellI]*TSSd_[cellI];
                }
            }
        }
    }
    Kss_.correctBoundaryConditions();
    S_.correctBoundaryConditions();
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HNGD::HNGD
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    hydrogenTransportSolver(mesh, mat, elementTransportDict, solverName),
    KD_
    (
        IOobject
        (
            "KD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    KN_
    (
        IOobject
        (
            "KN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    KG_
    (
        IOobject
        (
            "KG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Kss_
    (   
        IOobject
        (
            "Kss",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    S_
    (   
        IOobject
        (
            "S",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    hydrideChange_
    (   
        IOobject
        (
            "hydrideChange",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    t0_
    (
        IOobject
        (
            "t0",
            mesh_.time().timeName(),/*  */
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTime, -1),
        "calculated"
    ),
    V_H_("V_H", dimVolume/dimMoles, 1.7e-6) // (m3/mol)
{
    // Read options from dedicated subdict 
    dictionary HDict
    (
        elementTransportDict.
        subOrEmptyDict("HSolverOptions")
    );

    stressDrivenDiffusion_ = HDict.lookupOrDefault<bool>("stressDrivenDiffusion", false);

    //V_H_ = HDict.lookupOrDefault<dimensionedScalar>("V_H", 1.7e-6); // (m3/mol)
    if (HDict.found("V_H"))
    {
    #ifdef OPENFOAMFOUNDATION
        V_H_.value() = HDict.lookup<scalar>("V_H");
    #elif OPENFOAMESI
        V_H_.value() = HDict.get<scalar>("V_H");
    #endif
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void HNGD::correct()
{

    const dimensionedScalar& R = constant::physicoChemical::R;
    
    // Update all submodels
    hydrogenTransportSolver::correct();
    computeCtot();
    updateCoefficients();
    updateRates();

    const dictionary& stressControl =  mesh_.solutionDict().subDict("stressAnalysis");
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    int nInnerIter = 0;

    bool convergedInnerCss(false);
    bool convergedInnerCpp(false);

    // Interpolate gradT at faces
    surfaceVectorField gradTf = fvc::interpolate(gradT_);

    //- Build Soret term ----- [m2/s]
    surfaceScalarField phi
    (
        "phiH", 
        fvc::interpolate
        (
            diffH_ *  Q_ /(R*pow(T_, 2.0))
        ) * fvc::snGrad(T_) * mesh_.magSf()
    );

    // Stress driven diffusion term
    if(stressDrivenDiffusion_)
    {
        const volScalarField& sigmaHyd_ = mesh_.lookupObject<volScalarField>("sigmaHyd");
                
        phi -= fvc::interpolate
        (
           diffH_ * V_H_ /(R*T_)
        ) * fvc::snGrad(sigmaHyd_) * mesh_.magSf();
    }

    scalar residualCss = 0;
    scalar residualCpp = 0;
    scalar initialResidualCss;
    scalar initialResidualCpp;
    
    Css_.correctBoundaryConditions();
    
    volScalarField alpha
    (
        IOobject
        (
            "Kss_::relax",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("", dimless, 1.0)
    );
    
    do
    {
        // Store prev iter values for rel residual
        Css_.storePrevIter();
        Cpp_.storePrevIter();

        {
            // Equation for Concentration in Solid Solution
            fvScalarMatrix CssEqn
            (
                fvm::ddt(Css_)
              - fvm::laplacian(diffH_, Css_)
              - fvm::div(phi, Css_, "div(phi,Css)") // phi has contribution from the Soret term as well as the Stress driven diffusion(if present)
              - fvm::Sp(alpha*Kss_, Css_)
              - S_
            );

            initialResidual_ = CssEqn.solve().initialResidual();

            if (debug > 1)
            {
                surfaceScalarField flux = CssEqn.flux();

                Info<< "Hydrogen flux summary:" << nl;

                forAll(flux.boundaryField(), patchI)
                {
                    const scalarField& pf = flux.boundaryField()[patchI];

                    Info<< tab << mesh_.boundaryMesh()[patchI].name() << ": " << gSum(pf) << nl;
                }

                Info<< endl;
            }

            // Correct concentration in precipitation hydrides
            residualCpp = solve
            (
                fvm::ddt(Cpp_)
              + alpha*Kss_*Css_
              + S_
            ).initialResidual();
                        
            // Correct for overshoots
            const volScalarField& Css0 = Css_.oldTime();
            
            forAll(Css_, cellI)
            {
                scalar dC;
                
                if (Css0[cellI] > TSSd_[cellI])
                {
                    if (Css0[cellI] > Css_[cellI])
                    {
                        // Overshoot Tssd from growth region
                        alpha[cellI] = min((Css0[cellI] - TSSd_[cellI])
                                            / max(Css0[cellI] - Css_[cellI], SMALL)*alpha[cellI],
                                            1.0);
                        
                        dC = max(TSSd_[cellI] - Css_[cellI], 0);

                        Css_[cellI] += dC;
                        Cpp_[cellI] -= dC;
                    }
                    else
                    {
                        alpha[cellI] = 1;
                    }
                }
                else
                {
                    if (Css0[cellI] > Css_[cellI])
                    {
                        // Overshoot Tssd from dissolution region
                        alpha[cellI] = min((TSSd_[cellI] - Css0[cellI])
                                            / max((Css_[cellI] - Css0[cellI]), SMALL)*alpha[cellI],
                                            1.0);
                                           
                        dC = max(Css_[cellI] - TSSd_[cellI], 0);

                        Css_[cellI] -= dC;
                        Cpp_[cellI] += dC;
                    }
                    else
                    {
                        alpha[cellI] = 1;
                    }
                }
                
                //- Negative concentrations
                dC = (min(Css_[cellI], 0))
                          - (min(Cpp_[cellI], 0));

                Css_[cellI] -= dC;
                Cpp_[cellI] += dC;
            }
                        
            //hydrideChange needed for hydride reorientation model
            hydrideChange_ = Cpp_ - Cpp_.oldTime(); //Kss_*Css_ + S_;

            Css_.correctBoundaryConditions();
            Cpp_.correctBoundaryConditions();
            hydrideChange_.correctBoundaryConditions();
        }

        computeCtot();
        updateCoefficients();
        updateRates();

        // Store initial residual if in first inner iteration
        if(nInnerIter == 0)
        {
            initialResidualCss = residualCss; 
            initialResidualCpp = residualCpp; 
        }

        // Calculate a different residual based on the relative change -- Css
        scalar denomCss = gMax(mag(Css_.primitiveField() 
                                   - Css_.oldTime().primitiveField()));
        
        scalar numCss = gMax(mag(Css_.primitiveField() 
                                 - Css_.prevIter().internalField()));

        if (denomCss < SMALL)
        {
            denomCss = max(gMax(Css_.primitiveField()), SMALL);
        }
        
        scalar relResidualCss = numCss/denomCss;

        if(denomCss < 1e-5)
        {
            relResidualCss = VSMALL;
        }

        Info << "relResidualCss " << relResidualCss << endl;
        convergedInnerCss = converged(Css_.name());

        // Calculate a different residual based on the relative change -- Cpp
        scalar denomCpp = gMax(mag(Cpp_.primitiveField() 
                                   - Cpp_.oldTime().primitiveField()));
        
        scalar numCpp = gMax(mag(Cpp_.primitiveField() 
                                 - Cpp_.prevIter().internalField()));

        if (denomCpp < SMALL)
        {
             denomCpp = max(gMax(Cpp_.primitiveField()), SMALL);
        }
        
        scalar relResidualCpp_ = numCpp/max(denomCpp, SMALL);

        Info << "relResidualCpp " << relResidualCpp_ << endl;
        
    } while 
    (
        (not(convergedInnerCss)) 
        && (++nInnerIter < nCorr)
    );

    initialResidual_ = max(initialResidualCss, initialResidualCpp);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
