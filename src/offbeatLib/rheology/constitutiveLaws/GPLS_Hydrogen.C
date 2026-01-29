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

#include "GPLS_Hydrogen.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
#include "Time.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hydrogenGPLS, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        hydrogenGPLS, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
/////parameters hydrogenGPLS
Foam::scalar Foam::hydrogenGPLS::m
(
    const scalar& T
)
{
    const scalar mT = 1 / ( 1 + exp(10.2 * (T/292 -1)));
	const scalar m = 1 / (77.68 * mT + 4.11 * (1 - mT));
    return  m;
}

//- Strength modulus
Foam::scalar Foam::hydrogenGPLS::K
(
    const scalar& T,
	const scalar& Css,
    const scalar& Cpp
)
{
    //-- KT
    const scalar KT = 1 / (1 + exp(1.77 * (T/1007 - 1)));
	
	//-- K
	
    return (1 - 1.175e-4*Css + (6.15e-5 - 4.38e-8*T)*Cpp) *
			((1.409e9 - 8.952e5*T)*KT + 4.05e7 * (1 - KT));
}

Foam::scalar Foam::hydrogenGPLS::n0
(
    const scalar& T,
    const scalar& Cpp
)
{   
	//- n0T
	const scalar N0T = 1 / (1 + exp(12*(T/(810 - 9.19e-2*Cpp)-1))); 
    //-- n0
    return( (1 + 1.45e-4*Cpp) * (4.86e-2*N0T + 2.35e-2 * (1-N0T)) );
}

Foam::scalar Foam::hydrogenGPLS::alpha_n
(
    const scalar& T,
	const scalar& Cpp
)
{
    return (53.16 + 1.27e-2*Cpp) * (1 + exp(11.1*(T/738 - 1)));
}

Foam::scalar Foam::hydrogenGPLS::L
(
    const scalar& T,
    const scalar& Cpp,
    const scalar& p
)
{

    return pow(p+1e-4, n0(T, Cpp)) * exp( -alpha_n(T, Cpp) * p ) + ( 1 - exp(-alpha_n(T, Cpp) * p ));
}

// anistropy coefficients and Hill sigma equivalent
Foam::scalar Foam::hydrogenGPLS::Hrr
(
    const scalar& T
)
{
    return 0.485 + 9.5e-2 / ( 1 + exp( 12 * (T/740 - 1) ));
}

Foam::scalar Foam::hydrogenGPLS::Hzz
(
    const scalar& T
)
{
    return 0.52 + (-0.23 + 4e-4*T) / (1 + exp( 15 * (T/550 - 1) ));
}

Foam::scalar Foam::hydrogenGPLS::Htt
(
    const scalar& T
)
{   
    return 1 - Hrr(T);
}

// sigma Equivalent function
Foam::scalar Foam::hydrogenGPLS::sigmaEq_func
(
    const symmTensor& sigmaEffCyl,
    const scalar& T

)
{
    const scalar Hshear = 1.5;

    const scalar sigmaEq = sqrt(
      Hrr(T) * pow(sigmaEffCyl.yy() - sigmaEffCyl.zz(),2 ) 
    + Htt(T) * pow(sigmaEffCyl.zz() - sigmaEffCyl.xx(), 2) 
    + Hzz(T) * pow(sigmaEffCyl.xx() - sigmaEffCyl.yy(), 2) 
    + 2 * Hshear * pow(sigmaEffCyl.xy(), 2) 
    + 2 * Hshear * pow(sigmaEffCyl.xz(), 2) 
    + 2 * Hshear * pow(sigmaEffCyl.yz(), 2));
    return sigmaEq;
}

// cartesian -> cylindrical coordinates
Foam::symmTensor Foam::hydrogenGPLS::cylSig
(
    const symmTensor& sigma,
    const vector& Ci
)
{
    symmTensor sigmaEffCyl = sigma *  scalar(0);

    tensor R = tensor::zero;

    const scalar  radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);

    R.xx() =  Ci[0]/radius;
    R.xy() =  Ci[1]/radius;
    R.yy() =  Ci[0]/radius;
    R.yx() = -Ci[1]/radius;
    R.zz() = 1;
    sigmaEffCyl = symm( R & sigma & R.T() );

    return sigmaEffCyl;
}

Foam::scalar Foam::hydrogenGPLS::G
(
    const scalar& sigmaEq,
    const scalar& T,
    const scalar& Css,
    const scalar& Cpp,
    const scalar& epsilonPEq
)
{
    return pow(sigmaEq / (L(T, Cpp, epsilonPEq) * K(T, Css, Cpp)), 1/m(T));
}

Foam::scalar Foam::hydrogenGPLS::F
(
    const scalar& sigmaEq,
    const scalar& T,
    const scalar& Css,
    const scalar& Cpp,
    const scalar& epsilonPEq0,
    const scalar& dEpsilonPEq,
    const scalar& deltaT,
    const scalar& theta
)
{
    return dEpsilonPEq - deltaT * G(sigmaEq, T, Css, Cpp, epsilonPEq0+theta*dEpsilonPEq);
}

//- Derivative of F
Foam::scalar Foam::hydrogenGPLS::dF
(   
    const scalar& sigmaEq,
    const scalar& T,
    const scalar& Css,
    const scalar& Cpp,
    const scalar& p,
    const scalar& dEpsilonPEq,
    const scalar& deltaT,
    const scalar& theta
)
{
   // dEpsilonPEq = max(dEpsilonPEq, 0);
   const scalar one_m = 1 / m(T);
   const scalar pDot = G(sigmaEq, T, Css, Cpp, p + theta*dEpsilonPEq);

   const scalar term1 = deltaT * one_m * pDot / L(T, Cpp, p + theta*dEpsilonPEq);

   const scalar Lp = p + theta* dEpsilonPEq + 1e-4;
   const scalar Ln0 = pow(Lp, n0(T, Cpp));
   const scalar Lalpha = exp(-alpha_n(T, Cpp) * (p + theta*dEpsilonPEq));

   const scalar dL =( Ln0*Lalpha * (n0(T,Cpp)/Lp - alpha_n(T, Cpp)) ) + (Lalpha*alpha_n(T,Cpp));

   return 1 + theta * term1 * dL ;
}

/*
------ Hill Tensor nota -----
- Ideally we should modifiy OpenFOAM to add a new class of field... (see FOAM extend with the symmTensor4thOrder class) tedious work and would necessitate to touch the macros functions handlings fields in OpenFOAM. 
- I looked at it. It looks scary.
- We dont really need a full 4th order tensor class with all the mathematical operations implemented... just the double dot product!
- We create 6 (shear component are equals to 1.5, see Hshear_ volScalarField) scalarfield corresponding to each component of the symmetric 4th order hill tensor
  "xxxx", "xxyy", "xxzz",
          "yyyy", "yyzz",
                  "zzzz",
                         "xyxy",
                                "yzyz",
                                       "zxzx" 

------ Double dot product of a symm4thOrderTensor and a symmTensor -----
t4th = symm 4th order tensor
st = symm tensor
(
    t4th.xxxx()*st.xx() + t4th.xxyy()*st.yy() + t4th.xxzz()*st.zz(),
    t4th.xyxy()*st.xy(),
    t4th.zxzx()*st.xz(),

    t4th.xxyy()*st.xx() + t4th.yyyy()*st.yy() + t4th.yyzz()*st.zz(),
    t4th.yzyz()*st.yz(),
    t4th.xxzz()*st.xx() + t4th.yyzz()*st.yy() + t4th.zzzz()*st.zz()
);

- We dont create the 4th order tensor but we calculate each of its components
- and we hard code the (contracted) double dot product H : (sigmaEff/sigmaEq);
*/
Foam::symmTensor Foam::hydrogenGPLS::normal_tensor
(
    const scalar& T,
    const scalar& sigmaEq,
    const symmTensor& sigmaEff

)
{
    //- SigmaEff/sigmaEq
    const symmTensor sig_ratio = sigmaEff / sigmaEq;

    // - Scalar parameters Hill tensor (Lemaitre - Caboche formalism)
    const scalar F = Hzz(T);
    const scalar G = Hrr(T);
    const scalar H = Htt(T);
    
    const scalar Hxxxx = F + H;
    const scalar Hxxyy = -F;
    const scalar Hxxzz = -H;
    const scalar Hyyyy = G + F;
    const scalar Hyyzz = - G;
    const scalar Hzzzz = H + G;
    const scalar Hshear = 1.5;

    return  symmTensor(
            Hxxxx  * sig_ratio.xx() + Hxxyy*sig_ratio.yy() + Hxxzz*sig_ratio.zz(),

            Hshear * sqrt(2.0) * sig_ratio.xy(),
            Hshear * sqrt(2.0) *  sig_ratio.xz(), 

            Hxxyy  * sig_ratio.xx() + Hyyyy*sig_ratio.yy() + Hyyzz*sig_ratio.zz(),

            Hshear * sqrt(2.0)*sig_ratio.yz(),

            Hxxzz  * sig_ratio.xx() + Hyyzz*sig_ratio.yy() + Hzzzz*sig_ratio.zz()
        );
}

void Foam::hydrogenGPLS::correctHydrogenConcentrations
(
    const labelList& addr
)
{
    const volScalarField& T_(mesh_.lookupObject<volScalarField>("T"));
    const scalarField& Ti = T_.internalField();

    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];
        scalar TSS = 99000 * exp(-34523/8.314/Ti[cellI]);
        
        Css_[cellI] = min(Ctot_[cellI], TSS);
        Cpp_[cellI] = Ctot_[cellI] - Css_[cellI];

        const cell& c = mesh_.cells()[cellI];  
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            // forAll(patchIDs, registeredID)
            {
                if (patchID > -1 and Css_.boundaryField()[patchID].size())
                {

                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);           
        
                    const scalarField TP(T_.boundaryField()[patchID]);
                    scalarField& CssP(Css_.boundaryFieldRef()[patchID]);
                    scalarField& CppP(Css_.boundaryFieldRef()[patchID]);
                    scalarField& CtotP(Ctot_.boundaryFieldRef()[patchID]);
                    
                    scalar TSSP = 99000 * exp(-34523/8.314/TP[faceID]);
                    CssP[faceID] = min(CtotP[faceID], TSSP);
                    CppP[faceID] = CtotP[faceID] - CssP[faceID];
                }
            }
        }
    }
    //Css_.correctBoundaryConditions();
    //Cpp_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::hydrogenGPLS::hydrogenGPLS
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 0.1 )),
    theta_(lawDict.lookupOrDefault<scalar>("theta", 0.5)),
    tol_(lawDict.lookupOrDefault<scalar>("error", 1e-6)),
    maxIter_(lawDict.lookupOrDefault<scalar>("max_iter", 1500)),
    epsilonP_(createOrLookup<symmTensor>(mesh, "epsilonP")),
    DEpsilonP_(createOrLookup<symmTensor>(mesh, "DEpsilonP")),
    epsilonPEq_(createOrLookup<scalar>(mesh, "epsilonPEq")),
    DEpsilonPEq_(createOrLookup<scalar>(mesh, "DEpsilonPEq")),
    activeYield_(createOrLookup<scalar>(mesh, "activeYield")),
    plasticN_(createOrLookup<symmTensor>(mesh, "plasticN")),
    Css_(mesh.lookupObjectRef<volScalarField>("Css")),
    Cpp_(mesh.lookupObjectRef<volScalarField>("Cpp")),
    Ctot_(mesh.lookupObjectRef<volScalarField>("Ctot")),
    m_
    (
        IOobject
        (
            "m",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("m", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    Kmod_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("K", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    L_
    (
        IOobject
        (
            "L",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("L", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    etaVis_
    (
        IOobject
        (
            "etaVis",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("etaVis", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    sigmaEqH_
    (
        IOobject
        (
            "sigmaEqH",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sigmaEqH", dimPressure, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    maxAverageEpsilonPEq_(GREAT),
    maxMaximumEpsilonPEq_(GREAT)
{  
    const dictionary& controlDict(mesh_.time().controlDict());

    bool adjustTime = false;
    if (controlDict.found("adjustableTimeStep"))
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustableTimeStep", false);
    }
    else
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustTimeStep", false);
    }

    if(adjustTime)
    {
        controlDict.lookup("maxAverageEpsilonPEq") >> maxAverageEpsilonPEq_;
        controlDict.lookup("maxMaximumEpsilonPEq") >> maxMaximumEpsilonPEq_;
    }

    epsilonP_.oldTime();
    DEpsilonP_.oldTime();
    epsilonPEq_.oldTime();
    DEpsilonPEq_.oldTime();
    activeYield_.oldTime();
    plasticN_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hydrogenGPLS::~hydrogenGPLS()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hydrogenGPLS::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
	correctHydrogenConcentrations(addr);

    scalar i = 0;
    scalar error = 10;

    scalar dF_n = 0;
    scalar F_n = 0;

    //- TODO: Instead of looping over cells one might store the addressing and
    //- work with openfoam field syntax.
    const volScalarField& mu_(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& K_(mesh_.lookupObject<volScalarField>("threeK"));
    const volScalarField& lambda_(mesh_.lookupObject<volScalarField>("lambda"));
    const scalarField& muI = mu_.internalField();
    const scalarField& lambdaI = lambda_.internalField();

    //-Hydrogen concentration in solid solution and in precipitates fields
    const volScalarField& Css_(mesh_.lookupObject<volScalarField>("Css"));
    const volScalarField& Cpp_(mesh_.lookupObject<volScalarField>("Cpp"));

    // Take references to the internal fields for efficiency
    symmTensorField& plasticNI = plasticN_.ref();

    //- Take references to internal fields for fields that must be changed
    scalarField& DEpsilonPEqI = DEpsilonPEq_.ref();
    symmTensorField& DEpsilonPI = DEpsilonP_.ref();
    symmTensorField& epsilonI = epsilonEl.ref();

    symmTensorField& epsilonPI = epsilonP_.ref();
    scalarField& epsilonPEqI =  epsilonPEq_.ref();

    const symmTensorField& epsilonPOldTimeI = epsilonP_.oldTime().internalField(); 
    const scalarField& epsilonPEqOldTimeI = epsilonPEq_.oldTime().internalField();
    const scalarField& DEpsilonPEqOldTimeI = DEpsilonPEq_.oldTime().internalField();

    const scalar& deltaT = mesh_.time().deltaT().value();
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    //- Check if the simulation is large strain
    bool largeStrain = mesh_.foundObject<volTensorField>("F") ?
    true : false;

    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];

        //- Cell temperature
        const scalar TOld = Temp.oldTime().internalField()[cellI];
        const scalar T = Temp.internalField()[cellI]; 
        const scalar Tm = TOld + theta_*(T-TOld);

        const scalar CssOld = Css_.oldTime().internalField()[cellI];
        const scalar CssI = Css_.internalField()[cellI]; 
        const scalar Cssm = CssOld + theta_*(CssI-CssOld);

        const scalar CppOld = Cpp_.oldTime().internalField()[cellI];
        const scalar CppI = Cpp_.internalField()[cellI]; 
        const scalar Cppm = CppOld + theta_*(CppI-CppOld);

        // hydrogenGPLS PARAMETERS. For tracking and validation
        m_[cellI] = m(T);
        Kmod_[cellI] = K(T,  CssI, CppI);

        // Calculate trial stress
        const symmTensor sigmaEff = 2.0*muI[cellI]*epsilonI[cellI]
        + lambdaI[cellI]*tr(epsilonI[cellI])*I;

        //convert to cylindrical coordinates
        const vector& Ci = mesh_.C()[cellI];
        const symmTensor sigmaEffCyl = cylSig(sigmaEff, Ci);

        // Calculate Sigma Hill Equivalent
        sigmaEqH_[cellI] = sigmaEq_func(sigmaEffCyl, Tm);

        //Init plastic vector at 0
        plasticNI[cellI] = symmTensor::zero;

        //Init the equivalent strain increment over the current time step
        scalar DEpsilonPEq_n = 0;
        DEpsilonPEqI[cellI] = 0.0;

        //- Flow ---> Newton-Raphson iteration
        if(sigmaEqH_[cellI] > SMALL)
        {
            // Calculate normal vector / plastic multiplier N = H : (sigma/sigmaEq)
            plasticNI += normal_tensor(Tm, sigmaEqH_[cellI], sigmaEffCyl);

            // Initialize Newton-Raphson
            // If we assume that the time step is small enough 
            // ---> DEpsilonPEq(t) ~ DEpsilonPEq(t+dt) and DEpsilonPEq(t) is a good first estimation for the Newton-Raphson method
            DEpsilonPEq_n = DEpsilonPEqI[cellI];

            i = 0;
            error = 10;
            // Newton Raphson iterations

            do{    
                //F
                F_n = F(sigmaEqH_[cellI], Tm, Cssm, Cppm,
                               epsilonPEqOldTimeI[cellI], DEpsilonPEq_n, deltaT, theta_);
                // Info << "F_n --> " << F_n << endl;

                //dF
                dF_n = dF(sigmaEqH_[cellI], Tm, Cssm, Cppm,
                                 epsilonPEqOldTimeI[cellI], DEpsilonPEq_n, deltaT, theta_);

                //Newton Raphson estimation
                // np1 --> n + 1 with n the iteration of the newton rahpson method
                scalar DEpsilonPEq_np1 = DEpsilonPEq_n - F_n/dF_n;

                error = mag(DEpsilonPEq_np1 - DEpsilonPEq_n);
				
                //Update NR estimation with a relaxat ion factor
                DEpsilonPEq_n = relax_ * DEpsilonPEq_np1 + (1-relax_) * DEpsilonPEq_n;
                i++;

            }while(error > tol_ && i < maxIter_);

            if(i >= maxIter_)
            {
                Info << "WARNING GPLS INTEGRATOR" << endl;
                Info << "max iterations for the NR method reached (cell computation)" << endl; 
                Info << "With a value error of: " << error << endl;
                Info << "Number of iterations: " << i << endl;
                Info << "Value of the increment: " << DEpsilonPEq_n << endl;
                Info << "" << endl;
            }

            //Update DEpsilonPEq with value obtained with the Newton-Rahpson method
            DEpsilonPEqI[cellI] = DEpsilonPEq_n;
        }

        // Update hardening and viscosity with the new epsilon viscoplastic equivalent strain
        L_[cellI] = L(T, Cpp_[cellI], epsilonPEqOldTimeI[cellI]+DEpsilonPEqI[cellI]);
        etaVis_[cellI] =  Kmod_[cellI] * L_[cellI];

        // Update epsilonPEq
        epsilonPEqI[cellI] = epsilonPEqOldTimeI[cellI] + DEpsilonPEqI[cellI];

        // Update DEpsilonP
        DEpsilonPI[cellI] = plasticNI[cellI] * DEpsilonPEqI[cellI];

        // Update total plastic strain
        epsilonPI[cellI] = epsilonPOldTimeI[cellI] + DEpsilonPI[cellI];

        // Update elastic strain
        epsilonI[cellI] -= DEpsilonPI[cellI];

        // Calculate deviatoric stress
        sigmaDev[cellI] = dev(sigmaEff) - 2*muI[cellI]*DEpsilonPI[cellI];

        // Calculate the hydrostatic pressure
        const scalar trEpsilon = tr(epsilonI[cellI]);

        sigmaHyd[cellI] = 1.0/3.0*K_[cellI]*trEpsilon;

        const cell& c = mesh_.cells()[cellI];  
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            // forAll(patchIDs, registeredID)
            {
                if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                {

                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);           
    
                    // Take references to the boundary patch fields for efficiency
                    symmTensorField& sigmaDevP = 
                    sigmaDev.boundaryFieldRef()[patchID];

                    scalarField& sigmaHydP = 
                    sigmaHyd.boundaryFieldRef()[patchID];

                    symmTensorField& plasticNP = 
                    plasticN_.boundaryFieldRef()[patchID];

                    const scalarField& muP = 
                    mu_.boundaryField()[patchID];

                    const scalarField& lambdaP = 
                    lambda_.boundaryField()[patchID];

                    scalarField& DEpsilonPEqP = 
                    DEpsilonPEq_.boundaryFieldRef()[patchID];

                    symmTensorField& DEpsilonPP = 
                    DEpsilonP_.boundaryFieldRef()[patchID];

                    symmTensorField& epsilonElP = 
                    epsilonEl.boundaryFieldRef()[patchID];

                    symmTensorField& epsilonPP = 
                    epsilonP_.boundaryFieldRef()[patchID];

                    scalarField& epsilonPEqP = 
                    epsilonPEq_.boundaryFieldRef()[patchID];

                    const scalarField& KP = K_.boundaryField()[patchID];

                    const symmTensorField& epsilonPOldTimeP = 
                    epsilonP_.oldTime().boundaryField()[patchID]; 

                    const scalarField& epsilonPEqOldTimeP = 
                    epsilonPEq_.oldTime().boundaryField()[patchID]; 

                    const scalarField& DEpsilonPEqOldTimeP = 
                    DEpsilonPEq_.oldTime().boundaryField()[patchID];

                    const scalar& TP = Temp.boundaryField()[patchID][faceID];
                    const scalar& TPOld = Temp.oldTime().boundaryField()[patchID][faceID];
                    const scalar& TPm = TPOld + theta_*(TP - TPOld);

                    const scalar& CssP = Css_.boundaryField()[patchID][faceID];
                    const scalar& CssPOld = Css_.oldTime().boundaryField()[patchID][faceID];
                    const scalar& CssPm = CssPOld + theta_*(CssP - CssPOld);

                    const scalar& CppP = Cpp_.boundaryField()[patchID][faceID];
                    const scalar& CppPOld = Cpp_.oldTime().boundaryField()[patchID][faceID];
                    const scalar& CppPm = CppPOld + theta_*(CppP - CppPOld);

                    scalarField& sigmaEqHP = sigmaEqH_.boundaryFieldRef()[patchID];

                    scalarField& LP =
                    L_.boundaryFieldRef()[patchID];
					
					scalarField& KmodP =
                    Kmod_.boundaryFieldRef()[patchID];
                    
                    scalarField& etaVisP =
                    etaVis_.boundaryFieldRef()[patchID];

                    // Calculate trial stress with viscoplastic correction
                    const symmTensor sigmaEffP = 2.0*muP[faceID]*epsilonElP[faceID]
                    + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                    //convert to cylindrical coordinates
                    const vector& Ci = mesh_.C()[cellI];
                    const symmTensor sigmaEffCylP = cylSig(sigmaEffP, Ci);

                    // Hill sigma quivalent
                    sigmaEqHP[faceID] = sigmaEq_func(sigmaEffCylP, TPm);

                    //Init plastic vector at 0
                    plasticNP[faceID] = symmTensor::zero;

                    //Init the equivalent strain increment over the current time step
                    scalar DEpsilonPEq_n = 0;
                    DEpsilonPEqP[faceID] = 0.0;

                    scalar iP = 0;
                    scalar errorP = 10;
                    if(sigmaEqHP[faceID] > SMALL){

                        plasticNP[faceID] += normal_tensor(TPm, sigmaEqHP[faceID], sigmaEffCylP);

                        // Initialize Newton-Raphson
                        DEpsilonPEq_n = DEpsilonPEqP[faceID];

                        do{
                            //F
                            scalar F_n = F(sigmaEqHP[faceID], TPm, CssPm, CppPm,
                                           epsilonPEqOldTimeP[faceID], DEpsilonPEq_n, deltaT, theta_);

                            //dF
                            scalar dF_n = dF(sigmaEqHP[faceID], TPm, CssPm, CppPm,
                                             epsilonPEqOldTimeP[faceID], DEpsilonPEq_n, deltaT, theta_);

                            scalar DEpsilonPEq_np1 = DEpsilonPEq_n - F_n/dF_n;

                            errorP = mag(DEpsilonPEq_np1 - DEpsilonPEq_n);

                            DEpsilonPEq_n = (relax_) * DEpsilonPEq_np1 + (1-relax_) * DEpsilonPEq_n;
                            iP++;

                        }while(errorP > tol_ && iP < maxIter_);
                        if(iP >= maxIter_)
                        {
                            Info << "WARNING GPLS INTEGRATOR" << endl;
                            Info << "max iterations for the NR method reached (patch computation)" << endl; 
                            Info << "With a value error of: " << errorP << endl;
                            Info << "Number of iterations: " << iP << endl;
                            Info << "Value of the increment: " << DEpsilonPEq_n << endl;
                            Info << "" << endl;
                        }

                        DEpsilonPEqP[faceID] = DEpsilonPEq_n;

                    }

                    LP[faceID] = L(TP, CppP, 
                                 epsilonPEqOldTimeP[faceID]+DEpsilonPEqP[faceID]);
                    etaVisP[faceID] = KmodP[faceID] * LP[faceID];
                    
                     // Update equivalent total equivalent plastic strain
                    epsilonPEqP[faceID] = epsilonPEqOldTimeP[faceID] + DEpsilonPEqP[faceID];

                    // Update DEpsilonP
                    DEpsilonPP[faceID] = DEpsilonPEqP[faceID]*plasticNP[faceID];

                    // Update total plastic strain
                    epsilonPP[faceID] = 
                    epsilonPOldTimeP[faceID] + DEpsilonPP[faceID];

                    // Update elastic strain
                    epsilonElP[faceID] -= DEpsilonPP[faceID];

                    // Calculate deviatoric stress
                    sigmaDevP[faceID] = 
                    dev(sigmaEffP) - 2*muP[faceID]*DEpsilonPP[faceID];

                    // Calculate the hydrostatic pressure
                    const scalar trEpsilonP = tr(epsilonElP[faceID]);

                    sigmaHydP[faceID] = 1.0/3.0*KP[faceID]*trEpsilonP;       
                }

            }
        }
    }
}


void Foam::hydrogenGPLS::correctEpsilonEl
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];

        epsilonEl[cellI] -= epsilonP_.oldTime()[cellI];
        
        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            {
                if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);     

                    // Take references to the boundary patch fields        
                    symmTensorField& epsilonEl_p = 
                    epsilonEl.boundaryFieldRef()[patchID]; 
                    
                    const symmTensorField& epsilonPOldTime_p = 
                    epsilonP_.oldTime().boundaryField()[patchID]; 

                    epsilonEl_p[faceID] -= epsilonPOldTime_p[faceID];
                }
            }
        }
    }
};


void Foam::hydrogenGPLS::updateTotalFields
(
    const labelList& addr
)
{
    // Info<< nl << "Updating total accumulated fields" << endl;
    // Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_()) << endl;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

//        sigmaY_.ref()[cellI] += DSigmaY_.internalField()[cellI];

        if (DEpsilonPEq_.internalField()[cellI] > SMALL)
        {
            activeYield_.ref()[cellI] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.ref()[cellI] = 0.0;
        }

        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if 
            (
                patchID > -1 
                and DEpsilonPEq_.boundaryField()[patchID].size()
                and !activeYield_.boundaryField()[patchID].coupled()
            )
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                if (DEpsilonPEq_.boundaryField()[patchID][faceID] > SMALL)
                {
                    activeYield_.boundaryFieldRef()[patchID][faceID] = 1.0;
                }
                else
                {
                    activeYield_.boundaryFieldRef()[patchID][faceID] = 0.0;
                }

            }
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    activeYield_.correctBoundaryConditions();

    const int nTotalCells = returnReduce(mesh_.nCells(), sumOp<int>());

    Info<< "    " << numCellsYielding << " cells ("
        << 100.0*scalar(numCellsYielding)/scalar(nTotalCells)
        << "% of the cells in this material) are actively yielding"
        << nl << endl;
    
}


void Foam::hydrogenGPLS::correctAdditionalStrain
(
    volSymmTensorField& additionalStrain,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        additionalStrain[cellI] = epsilonP_.internalField()[cellI];
    }
}

Foam::scalar Foam::hydrogenGPLS::nextDeltaT
(
    const labelList& addr
)
{
    //- Store old value of deltaT
    scalar oldDeltaT = mesh_.time().deltaT().value();

    //- Prepare average and maximum DepsilonCreepEq from last time step
    scalar averageDEpsilonPEq(0);
    scalar maxDEpsilonPEq(0);
    scalar volTot(0);

    const scalarField& DEpsilonPEqI(DEpsilonPEq_());
    const scalarField& V(mesh_.V());

    //- Find maximum and average DEpsilonEq
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        averageDEpsilonPEq += DEpsilonPEqI[cellI]*V[cellI];
        volTot += V[cellI];

        maxDEpsilonPEq = max(maxDEpsilonPEq, DEpsilonPEqI[cellI]);
    }

    reduce(volTot, sumOp<scalar>());

    reduce(averageDEpsilonPEq, sumOp<scalar>());

    reduce(averageDEpsilonPEq, maxOp<scalar>());

    averageDEpsilonPEq /= volTot;

    //- Calculate average and maximum creep rate from last time step
    scalar averageEpsilonPEqRate(averageDEpsilonPEq/oldDeltaT);
    scalar maxEpsilonPEqRate(maxDEpsilonPEq/oldDeltaT);

    Info << "deltaT: averageDEpsilonPEq " << averageDEpsilonPEq << endl;
    Info << "deltaT: averageEpsilonPEq Rate " << averageEpsilonPEqRate << endl;
    Info << "deltaT: maxDEpsilonPEq " << maxDEpsilonPEq << endl;
    Info << "deltaT: maxEpsilonPEqRate " << maxEpsilonPEqRate << endl;

    return min
    (
        maxAverageEpsilonPEq_/max(averageEpsilonPEqRate, SMALL), 
        maxMaximumEpsilonPEq_/max(maxEpsilonPEqRate, SMALL)
    );


}

// ************************************************************************* //


