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

#include "HSolverEq.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HSolverEq, 0);
    addToRunTimeSelectionTable
    (
        HSolverEq, 
        HSolverEq, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Update material properties relative to the hydrogen transport 
void Foam::HSolverEq::updateCoefficients()
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
                const cell& c = mesh_.cells()[cellI];
                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    if (patchID > -1 and hydrideVolFrac_.boundaryField()[patchID].size())
                    {
                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                        scalarField& hydrideVolFracP(hydrideVolFrac_.boundaryFieldRef()[patchID]);
                        // scalarField& alphaP(alpha_.boundaryFieldRef()[patchID]);

                        const scalarField& Tp(T_.boundaryField()[patchID]);
                        const scalarField& CppP(Cpp_.boundaryFieldRef()[patchID]);
                        const scalarField& CtotP(Ctot_.boundaryFieldRef()[patchID]);
                        const scalarField& TSSdP(TSSd_.boundaryFieldRef()[patchID]);

                        hydrideVolFracP[faceID] = hydrideVolFrac(Tp[faceID], CppP[faceID], TSSdP[faceID]);
                    }
                }
                hydrideVolFrac_.correctBoundaryConditions();
            }

            //- Diffusion coefficient & TSSp & TSSd
            HDiffusionCoefficientModel_->correctDH(addr);
            TSSpModel_->correctTSSp(addr);
            TSSdModel_->correctTSSd(addr);
        }
    }
}

//- Rate terms that we add to the differential equations 
void Foam::HSolverEq::updateConcentration()
{
    const labelListList& matAddrList(mat_.matAddrList());
    forAll(mat_.materialsList(), i)
    {
         if(mat_.materialsList()[i].name()=="Zircaloy")
        {
            const labelList& addr(matAddrList[i]);
            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];

                if(Css_[cellI] > TSSp_[cellI] && hydrideVolFrac_[cellI] <= 0.99)
                {                   
                    Css_[cellI] = TSSp_[cellI];
                    Cpp_[cellI] = Ctot_[cellI] - TSSp_[cellI];
                }
                else if(TSSp_[cellI] >= Css_[cellI] && Css_[cellI] > TSSd_[cellI])
                {
                }
                //- Case 2 : Dissolution of Hydrides into solid solution 
                else if(Cpp_[cellI] > 0 && Css_[cellI] < TSSd_[cellI])
                { 
                    Css_[cellI] = TSSd_[cellI];
                    Cpp_[cellI] = Ctot_[cellI] - TSSd_[cellI];
                }
                else
                {
                }
            }
        }
    }
    Css_.correctBoundaryConditions();
    Cpp_.correctBoundaryConditions();
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*General Rq on the structure of the code.
For a reason I dont get, we have to solve for Cpp.
In theory, it should be possible to not solve for Cpp and to just update the concentrations
explicitely. I tried without sucess, the total hydrogen content when doing so is not conserved*/
Foam::HSolverEq::HSolverEq
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict
)
:
    hydrogenTransportSolver(mesh, mat, elementTransportDict),
    deltaConcentration_
    (
        IOobject
        (
            "deltaConcentration",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0),
        "calculated"
    )
{
    dictionary HDict
    (
        elementTransportDict.
        subOrEmptyDict("HSolverOptions")
    );        

    //- read kinetic flag
    //kineticParameter_ = HDict.lookupOrDefault<scalar>("kineticParameter", 0);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::HSolverEq::correct()
{
    // Update thermal diffusion coefficient
    updateCoefficients();

    const dictionary& stressControl =  mesh_.solutionDict().subDict("stressAnalysis");
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    int nInnerIter = 0;

    bool convergedInnerCss(false);
    bool convergedInnerCpp(false);

    Info << "total hydrogen content----->  " << gSum(Ctot_ * mesh_.V().field()) << endl;
    
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
        // Cpp_.storePrevIter()
        // Equation for Concentration in Solid Solution
        fvScalarMatrix CssEqn
        (
            fvm::ddt(Css_)
          - fvm::laplacian(diffH_,Css_)
          - fvc::div(soretTerm * mesh_.Sf() & gradTf )
        );  
        
        scalar residualCss_ = CssEqn.solve().max().initialResidual();
        CssEqn.relax();
        Css_.relax();

        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            const scalar& initialResidualCss_ = residualCss_; 
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

    } while 
    (
        not(convergedInnerCss) 
        && ++nInnerIter < nCorr 
    );
    computeCtot();
    updateConcentration();
}

// ************************************************************************* //
