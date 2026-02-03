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

#include "HSolverOneTSS.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HSolverOneTSS, 0);
    addToRunTimeSelectionTable
    (
        HSolverOneTSS, 
        HSolverOneTSS, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Update material properties relative to the hydrogen transport 
void Foam::HSolverOneTSS::updateCoefficients()
{
    const labelListList& matAddrList(mat_.matAddrList());
    forAll(mat_.materialsList(), i)
    {
        if(mat_.materialsList()[i].name()=="Zircaloy")
        {
            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                //- VOLUME FRACTION OF HYDRIDES
                hydrideVolFrac_[cellI] = hydrideVolFrac(T_[cellI], Cpp_[cellI], TSSd_[cellI]);

                //- PRECIPITATION KINETIC PARAMETER
                if(kineticParameter_ == 0)
                {
                    alpha_[cellI] = pow(62.3,2) * exp(-2*4.12e4 / (8.31446 * T_[cellI]));
                }
                else
                {
                    alpha_[cellI] = kineticParameter_;
                }

                const cell& c = mesh_.cells()[cellI];
                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    if (patchID > -1 and hydrideVolFrac_.boundaryField()[patchID].size())
                    {
                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                        scalarField& hydrideVolFracP(hydrideVolFrac_.boundaryFieldRef()[patchID]);
                        scalarField& alphaP(alpha_.boundaryFieldRef()[patchID]);

                        const scalarField& Tp(T_.boundaryField()[patchID]);
                        const scalarField& CppP(Cpp_.boundaryFieldRef()[patchID]);
                        const scalarField& CtotP(Ctot_.boundaryFieldRef()[patchID]);
                        const scalarField& TSSdP(TSSd_.boundaryFieldRef()[patchID]);



                        hydrideVolFracP[faceID] = hydrideVolFrac(Tp[faceID], CppP[faceID], TSSdP[faceID]);
                        if(kineticParameter_ == 0)
                        {
                            alphaP[faceID] = pow(62.3,2) * exp(-2*4.12e4 / (8.31446 * Tp[faceID])); 
                        }
                        else
                        {
                            alphaP[faceID] = kineticParameter_;
                        }
                    }
                }
                hydrideVolFrac_.correctBoundaryConditions();
                alpha_.correctBoundaryConditions();
            }

            //- Diffusion coefficient & TSSp & TSSd
            HDiffusionCoefficientModel_->correctDH(addr);
            //TSSpModel_->correctTSSp(addr);
            TSSdModel_->correctTSSd(addr);
        }
    }
}

//- Rate terms that we add to the differential equations 
void Foam::HSolverOneTSS::updateRateTerms()
{
    const labelListList& matAddrList(mat_.matAddrList());
    forAll(mat_.materialsList(), i)
    {
         if(mat_.materialsList()[i].name()=="Zircaloy")
        {
            const labelList& addr(matAddrList[i]);
            forAll(addr, j)
            {
                const label cellI = addr[j];
                CssEqnAdd_[cellI] = 0;
                CppEqnAdd_[cellI] = 0;

                //- Case 1 :precipitation of hydrogen ---- Css > TSSP
                if(Css_[cellI] > TSSd_[cellI] && hydrideVolFrac_[cellI] <= 0.99)
                {                   
                    CppEqnAdd_[cellI] = alpha_[cellI] * (Css_[cellI] - TSSd_[cellI]);
                    CssEqnAdd_[cellI] = - CppEqnAdd_[cellI];
                }
                //- Case 2 : Dissolution of Hydrides into solid solution 
                else if(Cpp_[cellI] > 0 && Css_[cellI] <= TSSd_[cellI])
                {
                    CssEqnAdd_[cellI] = min(alpha_[cellI]*100,1) * (TSSd_[cellI] - Css_[cellI]);
                    CppEqnAdd_[cellI] = - CssEqnAdd_[cellI];  
                }
                else
                {
                    CssEqnAdd_[cellI] = 0;
                    CppEqnAdd_[cellI] = 0;
                }

                const cell& c = mesh_.cells()[cellI];  
                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);                  
                    if (patchID > -1 and Css_.boundaryField()[patchID].size())
                    {
                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                        const scalarField& hydrideVolFracP(hydrideVolFrac_.boundaryFieldRef()[patchID]);

                        const scalarField& CssP(Css_.boundaryFieldRef()[patchID]);
                        const scalarField& CppP(Css_.boundaryFieldRef()[patchID]);
                        const scalarField& TSSdP(TSSd_.boundaryFieldRef()[patchID]);

                        const scalarField& alphaP(alpha_.boundaryFieldRef()[patchID]);

                        scalarField& CssEqnAddP(CssEqnAdd_.boundaryFieldRef()[patchID]);
                        scalarField& CppEqnAddP(CppEqnAdd_.boundaryFieldRef()[patchID]);

                        CssEqnAddP[faceID] = 0;
                        CppEqnAddP[faceID] = 0;

                        if(CssP[faceID] > TSSdP[faceID] && hydrideVolFracP[faceID] <= 0.99)
                        {   
                            CppEqnAddP[faceID] = alphaP[faceID] * (CssP[faceID] - TSSdP[faceID]);
                            CssEqnAddP[faceID] = -CppEqnAddP[faceID];
                        }
                        else if(CppP[faceID] > 0 && CssP[faceID] <= TSSdP[faceID])
                        {
                            CssEqnAddP[faceID] = min(alphaP[faceID]*100,1) * (TSSdP[faceID] - CssP[faceID]);
                            CppEqnAddP[faceID] = -CssEqnAddP[faceID];
                        }
                        else
                        {
                            CssEqnAddP[faceID] = 0;
                            CppEqnAddP[faceID] = 0;
                        }
                    }
                }
            }
        }
    }
    CssEqnAdd_.correctBoundaryConditions();
    CppEqnAdd_.correctBoundaryConditions();
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*General Rq on the structure of the code.
For a reason I dont get, we have to solve for Cpp.
In theory, it should be possible to not solve for Cpp and to just update the concentrations
explicitely. I tried without sucess, the total hydrogen content when doing so is not conserved*/
Foam::HSolverOneTSS::HSolverOneTSS
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict
)
:
    hydrogenTransportSolver(mesh, mat, elementTransportDict),
    alpha_
    (
        IOobject
        (
            "alphaDiff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimensionSet(0, 0, -0.5, 0, 0, 0, 0), 1.0),
        "calculated"
    )
{
    dictionary HDict
    (
        elementTransportDict.
        subOrEmptyDict("HSolverOneTSSOptions")
    );        

    //- read kinetic flag
    kineticParameter_ = HDict.lookupOrDefault<scalar>("kineticParameter", 0);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::HSolverOneTSS::correct()
{
    // Update thermal diffusion coefficient
    updateCoefficients();

    const dictionary& stressControl =  mesh_.solutionDict().subDict("stressAnalysis");
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    int nInnerIter = 0;

    bool convergedInnerCss(false);
    bool convergedInnerCpp(false);

    Info << "total hydrogen content----->  " << gSum(Cpp_ * mesh_.V().field()) + gSum(Css_*mesh_.V().field()) << endl;
    
    // Interpolate gradT at faces
    surfaceVectorField gradTf = fvc::interpolate(gradT_);

    //- Build Soret term ----- [m2/s]
    surfaceScalarField soretTerm = fvc::interpolate
    (
        diffH_ * Css_  *  Q_ /(R_*pow(T_, 2.0))
    );

    do{
        // Store prev iter values for rel residual
        Css_.storePrevIter();
        Cpp_.storePrevIter();

        // Equation for Concentration in Solid Solution
        fvScalarMatrix CssEqn
        (
            fvm::ddt(Css_)
          - fvm::laplacian(diffH_,Css_)
          - fvc::div(soretTerm * mesh_.Sf() & gradTf )
          - CssEqnAdd_
        );  
        
        scalar residualCss_ = CssEqn.solve().max().initialResidual();
        CssEqn.relax();
        Css_.relax();

        // Equation for Concentration in Hydrides Precipitates. Most stupidest equation ever 
        fvScalarMatrix CppEqn
        (
            fvm::ddt(Cpp_)
          - CppEqnAdd_
       );  

         // Compute residual
        scalar residualCpp_ = CppEqn.solve().max().initialResidual();
        CppEqn.relax();
        Cpp_.relax();

        computeCtot();
        updateRateTerms();

        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            const scalar& initialResidualCss_ = residualCss_; 
            const scalar& initialResidualCpp_ = residualCpp_; 
        }
        
        // Calculate a different residual based on the relative change -- Css
        scalar denomCss = 
        gMax(mag(Css_.primitiveField() - Css_.oldTime().primitiveField()));
        
        scalar numCss = 
        gMax(mag(Css_.primitiveField() - Css_.prevIter().internalField()));

        if (denomCss < SMALL)
        {
            denomCss = max(gMax(mag(Css_.primitiveField())), SMALL);
        }
        
        scalar relResidualCss_ = numCss/denomCss;

        if(denomCss < 1e-5)
        {
            relResidualCss_ = VSMALL;
        }

        Info << "relResidualCss " << relResidualCss_ << endl;
        convergedInnerCss = converged(Css_.name());

        // // Calculate a different residual based on the relative change -- Cpp
        scalar denomCpp = 
        gMax(mag(Cpp_.primitiveField() - Cpp_.oldTime().primitiveField()));
        
        scalar numCpp = 
        gMax(mag(Cpp_.primitiveField() - Cpp_.prevIter().internalField()));

        if (denomCpp < SMALL)
        {
             denomCpp = max(gMax(mag(Cpp_.primitiveField())), SMALL);
        }
        
         scalar relResidualCpp_ = numCpp/denomCpp;

        if(denomCpp < 1e-5)
        {
             relResidualCpp_ = VSMALL;
        }

        Info << "relResidualCpp " << relResidualCpp_ << endl;
        convergedInnerCpp = converged(Cpp_.name());
    
    } while 
    (
        (not(convergedInnerCss) || not(convergedInnerCpp))
        && ++nInnerIter < nCorr 
    );
}

// ************************************************************************* //
