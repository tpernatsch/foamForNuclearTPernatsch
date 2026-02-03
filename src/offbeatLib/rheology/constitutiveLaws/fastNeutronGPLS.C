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

#include "fastNeutronGPLS.H"
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
    defineTypeNameAndDebug(fastNeutronGPLS, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        fastNeutronGPLS, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
/////parameters fastNeutronGPLS
Foam::scalar Foam::fastNeutronGPLS::phiPhi
(
    const scalar& phi
)const
{
    scalar phiPhi = 1 - exp(-0.3 * phi);
    return phiPhi;
}

Foam::scalar Foam::fastNeutronGPLS::m
(
    const scalar& T,
    const scalar& phi
)const
{
    const scalar m0T = 1 / (1 + exp(10.2 * (T / 692 - 1)));
    const scalar m0 = 77.68 * m0T + (4.11 * (1 - m0T));
    const scalar mPhi = 0.46 * phiPhi(phi) / (1 + exp(14.49 * (T / 870 - 1)));
    const scalar m = 1 / (m0 * (1 - mPhi));
    return  m;
}

Foam::scalar Foam::fastNeutronGPLS::K
(
    const scalar& T,
    const scalar& phi,
    const scalar& zeta
)
{
    //-- K0T
    const scalar K0T = 1 / (1 + exp(6.62 * (T/1007 - 1)));
    //-- K0
    const scalar term1 = 1.409e9 - 1.1558e6*T;
    const scalar K0 = term1*K0T + (4.05e7 * (1 - K0T));
    //-- KPhi
    const scalar KPhi = (0.42 * (1 - 6.44e-2 * phi * zeta) 
    * phiPhi(phi)) / (1 + exp(19 * (T/768 - 1)));
    //-- K
    return (K0 * (1 + KPhi));

}

Foam::scalar Foam::fastNeutronGPLS::n0
(
    const scalar& T,
    const scalar& phi
)
{   
    //-- N0
    const scalar N0 = 1 / (1 + exp(4 * (T/620 - 1))); 
    //n0
    return (7.9e-2*N0 + 3e-2 * (1-N0)) 
    * (1 - ((0.24*phiPhi(phi)) / (1 + exp(7 * (T/630 - 1)))));
}

Foam::scalar Foam::fastNeutronGPLS::ninf
(
    const scalar& T,
    const scalar& phi
)
{
    //-- Ninf
    const scalar Ninf = 1 / (1 + exp(5 * (T/590  - 1)));   
    //-- ninf
    return(5.1e-2*Ninf + 5.6e-3 * (1-Ninf)) 
    * (1 - (0.63*phiPhi(phi)) / (1 + exp(7 * (T / 630 - 1))));     
}

Foam::scalar Foam::fastNeutronGPLS::alpha_n
(
    const scalar& T
)
{
    return 40.45 * exp(2.03e-3 * T);   
}

Foam::scalar Foam::fastNeutronGPLS::L
(
    const scalar& T,
    const scalar& phi,
    const scalar& p
)
{
    //-- L
    // return pow(p + 1e-4, n0(T, phi)) * exp(-alpha_n(T) * p) 
    //+ (pow(p + 1e-4, ninf(T, phi)) * (1 - exp(-alpha_n(T) * p)));
    return pow(mag(p) + 1e-4, n0(T, phi)) * exp(-alpha_n(T) * mag(p)) 
    + (pow(mag(p) + 1e-4, ninf(T, phi)) * (1 - exp(-alpha_n(T) * mag(p))));
}

// anistropy coefficients
Foam::scalar Foam::fastNeutronGPLS::Hrr
(
    const scalar& T,
    const scalar& phi
)
{
    return 0.485 + 9.5e-2 * (1- phiPhi(phi)) / (1 + exp(12 * (T/740 -1))) 
    + 0.32 * (phiPhi(phi) / (1 + exp(10 * (T/760 -1))));
}

Foam::scalar Foam::fastNeutronGPLS::Hzz
(
    const scalar& T,
    const scalar& phi
)
{
    return 0.52 + (-0.23 + 4e-4*T) * (1- phiPhi(phi)) 
    / (1 + exp(15 * (T/550 -1))) - 0.16 * (phiPhi(phi) / (1 + exp(20 * (T/920 -1))));
}

Foam::scalar Foam::fastNeutronGPLS::Htt
(
    const scalar& T,
    const scalar& phi
)
{   
    return 1 - Hrr(T, phi);
}

// sigma Equivalent function
Foam::scalar Foam::fastNeutronGPLS::sigmaEq_func
(
    const symmTensor& sigmaEffCyl,
    const scalar& T,
    const scalar& phi

)
{
    const scalar Hshear = 1.5;

    const scalar sigmaEq = sqrt(
      Hrr(T, phi) * pow(sigmaEffCyl.yy() - sigmaEffCyl.zz(),2 ) 
    + Htt(T, phi) * pow(sigmaEffCyl.zz() - sigmaEffCyl.xx(), 2) 
    + Hzz(T, phi) * pow(sigmaEffCyl.xx() - sigmaEffCyl.yy(), 2) 
    + 2 * Hshear * pow(sigmaEffCyl.xy(), 2) 
    + 2 * Hshear * pow(sigmaEffCyl.xz(), 2) 
    + 2 * Hshear * pow(sigmaEffCyl.yz(), 2));
    return sigmaEq;
}

// cartesian -> cylindrical coordinates
Foam::symmTensor Foam::fastNeutronGPLS::cylSig
(
    const symmTensor& sigma,
    const vector& Ci
)
{
    symmTensor sigmaEffCyl = sigma *  scalar(0);

    // Transform the sigmaEffCyl tensor in cylindrical coordinates
    tensor R = tensor::zero;

    const scalar  radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);

    R.xx() =  Ci[0]/radius;
    R.xy() =  Ci[1]/radius;
    R.yy() =  Ci[0]/radius;
    R.yx() = -Ci[1]/radius;
    R.zz() = 1;
    sigmaEffCyl = symm( R & sigma & R.T() );
    
    const scalar theta = atan(Ci[1]/Ci[0]);
    return sigmaEffCyl;
}

// viscoplastic strain rate function of the fastNeutronGPLS
Foam::scalar Foam::fastNeutronGPLS::G
(
    const scalar& sigmaEq,
    const scalar& T,
    const scalar& phi,
    const scalar& zeta,
    const scalar& epsilonPEq
)
{
    return pow(sigmaEq/(L(T, phi, epsilonPEq)*K(T, phi, zeta)), 1/m(T, phi));
}

// F function for the NR iterations. We search the 0 of this function!
Foam::scalar Foam::fastNeutronGPLS::F
(
    const scalar& sigmaEq,
    const scalar& T,
    const scalar& phi,
    const scalar& zeta,
    const scalar& epsilonPEq0,
    const scalar& dEpsilonPEq,
    const scalar& deltaT,
    const scalar& theta
)
{
    return dEpsilonPEq - deltaT * G(sigmaEq, T, phi, zeta, epsilonPEq0+theta*dEpsilonPEq);
}

// Derivative of F
Foam::scalar Foam::fastNeutronGPLS::dF
(   
    const scalar& sigmaEq,
    const scalar& T,
    const scalar& phi,
    const scalar& zeta,
    const scalar& p,
    const scalar& dEpsilonPEq,
    const scalar& deltaT,
    const scalar& theta
)
{
   const scalar one_m = 1 / m(T, phi);
   const scalar pDot = G(sigmaEq, T, phi, zeta, p + theta*dEpsilonPEq);

   const scalar term1 = deltaT * one_m * pDot / L(T, phi, p + theta*dEpsilonPEq);

   const scalar Lp = p + theta*dEpsilonPEq + 1e-4;
   const scalar Ln0 = pow(Lp, n0(T, phi));
   const scalar Lninf = pow(Lp, ninf(T, phi));
   const scalar Lalpha = exp(-alpha_n(T) * (p + theta*dEpsilonPEq));
    
   const scalar dL = Ln0 * Lalpha * (n0(T, phi) / Lp - alpha_n(T)) +
                     Lninf * (1-Lalpha) * (ninf(T, phi) / Lp - alpha_n(T)) +
                     alpha_n(T)*Lninf;

                     
   return 1 + term1 * dL * theta;
}

/*
------ Hill Tensor nota -----
- Ideally we should modifiy OpenFOAM to add a new class of field... 
  (see FOAM extend with the symmTensor4thOrder class) 
- We dont really need a full 4th order tensor class
- We just need the double dot product!
- We use 6 (shear component are equals to 1.5) scalars
- The symmetric 4th order Hill Tensor can be visualized as:
  "Hxxxx", "Hxxyy", "Hxxzz",
           "Hyyyy", "Hyyzz",
                    "Hzzzz",
                            "Hxyxy",
                                    "Hyzyz",
                                            "Hzxzx" 

Its components are expressed with 6 anisotropy coefficients: F, G, H, L, M, N.
L, M, N are the shear components that are here equal to 1.5 (isotropic case)
We have (Lemaître and Chaboche, Mechanics of solid materials, 1994):
Hxxxx = F + H      Hyyzz = - G
Hxxyy = - F        Hzzzz = H + G
Fxxzz = -H
Hyyyy = G + F

F, G, H are equal to Hzz, Hrr or Htt (model parameters)

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


*/
Foam::symmTensor Foam::fastNeutronGPLS::normal_tensor
(
    const scalar& T,
    const scalar& phi,
    const scalar& sigmaEq,
    const symmTensor& sigmaEff

)
{
    //- SigmaEff/sigmaEq
    const symmTensor sig_ratio = sigmaEff / sigmaEq;

    // - Scalar parameters Hill tensor (Lemaitre - Caboche formalism)
    const scalar F = Hzz(T, phi);
    const scalar G = Hrr(T, phi);
    const scalar H = Htt(T, phi);
    
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fastNeutronGPLS::fastNeutronGPLS
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict),
    creepModel_(creepModel::New(mesh, lawDict)),
    phiValue_(readScalar(lawDict.lookup("phi"))),
    zetaValue_(readScalar(lawDict.lookup("zeta"))),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 0.1 )),
    theta_(lawDict.lookupOrDefault<scalar>("theta", 0.5)),
    tol_(lawDict.lookupOrDefault<scalar>("error", 1e-8)),
    maxIter_(lawDict.lookupOrDefault<scalar>("max_iter", 2500)),
    epsilonP_(createOrLookup<symmTensor>(mesh, "epsilonP")),
    DEpsilonP_(createOrLookup<symmTensor>(mesh, "DEpsilonP")),
    epsilonPEq_(createOrLookup<scalar>(mesh, "epsilonPEq")),
    DEpsilonPEq_(createOrLookup<scalar>(mesh, "DEpsilonPEq")),
    plasticN_(createOrLookup<symmTensor>(mesh, "plasticN")),
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
        dimensionedScalar("m", dimless, 0)
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
        dimensionedScalar("K", dimless, 0)
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
        dimensionedScalar("L", dimless, 0)
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
        dimensionedScalar("etaVis", dimless, 0)
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
        dimensionedScalar("sigmaEqH", dimPressure, 0)
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
    plasticN_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fastNeutronGPLS::~fastNeutronGPLS()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fastNeutronGPLS::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    creepModel_->correctCreep(epsilonEl, addr);

    //- TODO: Instead of looping over cells one might store the addressing and
    //- work with openfoam field syntax.
    const volScalarField& mu_(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& K_(mesh_.lookupObject<volScalarField>("threeK"));
    const volScalarField& lambda_(mesh_.lookupObject<volScalarField>("lambda"));
    
    const scalarField& muI = mu_.internalField();
    const scalarField& lambdaI = lambda_.internalField();
    
    symmTensorField& plasticNI = plasticN_.ref();
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

        /// fastNeutronGPLS PARAMETERS. For tracking and validation
        m_[cellI] = m(T, phiValue_);
        Kmod_[cellI] = K(T, phiValue_, zetaValue_);
        //L at the "begining" of the time step is calculated based on the old value of the equivalent strain
        L_[cellI] = L(T, phiValue_, epsilonPEqOldTimeI[cellI]); 
        etaVis_[cellI] =  Kmod_[cellI] * L_[cellI];

        // Calculate trial stress
        const symmTensor sigmaEff = 2.0*muI[cellI]*epsilonI[cellI]
        + lambdaI[cellI]*tr(epsilonI[cellI])*I;

        //convert to cylindrical coordinates
        const vector& Ci = mesh_.C()[cellI];
        const symmTensor sigmaEffCyl = cylSig(sigmaEff, Ci);

        // Calculate Sigma Hill Equivalent
        sigmaEqH_[cellI] = sigmaEq_func(sigmaEffCyl, Tm, phiValue_);

        //Init plastic vector at 0
        plasticNI[cellI] = symmTensor::zero;

        //Init the equivalent strain increment over the current time step
        scalar DEpsilonPEq_n = 0;
        DEpsilonPEqI[cellI] = 0.0;

        //- Flow ---> Newton-Raphson iteration
        if(sigmaEqH_[cellI] > 1.){

            // Calculate normal vector / plastic multiplier N = H : (sigma/sigmaEq)
            plasticNI[cellI] += normal_tensor(Tm, phiValue_, sigmaEqH_[cellI], sigmaEff);

            // Initialize Newton-Raphson
            // If we assume that the time step is small enough 
            // ---> DEpsilonPEq(t) ~ DEpsilonPEq(t+dt) and DEpsilonPEq(t) is a good first estimation for the Newton-Raphson method
            //DEpsilonPdEqI[cellI] = DEpsilonPEqOldTimeI[cellI];
            DEpsilonPEq_n = DEpsilonPEqI[cellI];

            scalar i = 0;
            scalar error = 10;
            // Newton Raphson iterations
            do{    
                //F
                scalar F_n = F(sigmaEqH_[cellI], Tm, phiValue_, zetaValue_, 
                            epsilonPEqOldTimeI[cellI], DEpsilonPEq_n, deltaT, theta_);

                //dF
                scalar dF_n = dF(sigmaEqH_[cellI], Tm, phiValue_, zetaValue_, 
                            epsilonPEqOldTimeI[cellI], DEpsilonPEq_n, deltaT, theta_);

                //Newton Raphson estimation
                // np1 --> n + 1 with n the iteration of the newton rahpson method
                scalar DEpsilonPEq_np1 = DEpsilonPEq_n - F_n/dF_n;

                // Error d
                error = mag(DEpsilonPEq_np1 - DEpsilonPEq_n);

                //Update NR estimation with a relaxation factor
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
        L_[cellI] = L(T, phiValue_, epsilonPEqOldTimeI[cellI]+DEpsilonPEqI[cellI]);
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

                    scalarField& LP =
                    L_.boundaryFieldRef()[patchID];
                    
                    scalarField& KmodP =
                    Kmod_.boundaryFieldRef()[patchID];
                    
                    scalarField& etaVisP =
                    etaVis_.boundaryFieldRef()[patchID];

                    scalarField& sigmaEqHP = 
                    sigmaEqH_.boundaryFieldRef()[patchID];

                    // Calculate trial stress
                    const symmTensor sigmaEffP = 2.0*muP[faceID]*epsilonElP[faceID]
                    + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                    //convert to cylindrical coordinates
                    const vector& Ci = mesh_.boundaryMesh()[patchID].faceCentres()[faceID];
                    const symmTensor sigmaEffCylP = cylSig(sigmaEffP, Ci);

                    // Hill sigma quivalent
                    sigmaEqHP[faceID] = sigmaEq_func(sigmaEffCylP, TPm, phiValue_);

                    //Init plastic vector at 0
                    plasticNP[faceID] = symmTensor::zero;

                    //Init the equivalent strain increment over the current time step
                    scalar DEpsilonPEq_n = 0;
                    DEpsilonPEqP[faceID] = 0.0;

                    scalar iP = 0;
                    scalar errorP = 10;
                    if(sigmaEqHP[faceID] > SMALL){

                        plasticNP[faceID] += normal_tensor(TPm, phiValue_, sigmaEqHP[faceID], sigmaEffP);

                        // Initialize Newton-Raphson
                        DEpsilonPEq_n = DEpsilonPEqP[faceID];
                        do{
                            //F
                            scalar F_n = F(sigmaEqHP[faceID], TPm, phiValue_, zetaValue_, 
                            epsilonPEqOldTimeP[faceID], DEpsilonPEq_n, deltaT, theta_);

                            //dF
                            scalar dF_n = dF(sigmaEqHP[faceID], TPm, phiValue_, 
                            zetaValue_, epsilonPEqOldTimeP[faceID], DEpsilonPEq_n, deltaT, theta_);


                             //- New estimation
                            scalar DEpsilonPEq_np1 = DEpsilonPEq_n - F_n/dF_n;

                            errorP = mag(DEpsilonPEq_np1 - DEpsilonPEq_n);

                            //- relaxation
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

                    LP[faceID] = L(TP, phiValue_, epsilonPEqOldTimeP[faceID]+DEpsilonPEqP[faceID]);
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


void Foam::fastNeutronGPLS::correctEpsilonEl
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

void Foam::fastNeutronGPLS::correctAxialStrain
(
    volSymmTensorField& epsilon,
    const labelList& addr
)
{
    creepModel_->correctAxialStrain(epsilon, addr);
}

void Foam::fastNeutronGPLS::correctAdditionalStrain
(
    volSymmTensorField& additionalStrain,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        additionalStrain[cellI] = creepModel_->additionalStrain()[cellI]
                                + epsilonP_.internalField()[cellI];
    }
}

Foam::scalar Foam::fastNeutronGPLS::nextDeltaT
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


