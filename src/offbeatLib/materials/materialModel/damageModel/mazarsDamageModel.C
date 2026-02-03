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

#include "mazarsDamageModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
#include "Time.H"
#include "tensor2D.H"
#include "eig3.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mazarsDamageModel, 0);
    addToRunTimeSelectionTable
    (
        damageModel, 
        mazarsDamageModel, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::mazarsDamageModel::At
(
    const scalar& T
)
{
    return 3.5;
}


Foam::scalar Foam::mazarsDamageModel::Ac
(
    const scalar& T
)
{
    // Parameters by M. Reymond
    if(T<1373)
    {
        return 79.276 * pow(T, -0.611434);
    }
    else if(T>=1373 && T<1623)
    {
        return -0.004 * T + 6.492;
    }
    else
    {
        return 0;
    }

    // // Parameters by J.S.
    // const scalar a = 8.9719e-7;
    // const scalar b = -3.11509e-3;
    // const scalar c = 2.8357;
    // const scalar d = 9.47303e6;
    // const scalar e = -2.6154;
    // if (T < 1500)
    // {
    //     return a*T*T + b*T + c;
    // }
    // else
    // {
    //   return d*pow(T,e);
    // }
}


Foam::scalar Foam::mazarsDamageModel::Bt
(
    const scalar& T
)
{
    return 1750;
}


Foam::scalar Foam::mazarsDamageModel::Bc
(
    const scalar& T
)
{
    // Parameters by M. Reymond
    if(T < 1373)
    {
        return 8.9692e-5*pow(T,2) - 0.034*T + 181.10163;
    }
    else if(T>=1373 && T<1623)
    {
        return -0.8*T + 1398.4;
    }
    else
    {
        return  0.461538461538462 * T - 649.076923076923;
    }

    // // Parameters by J. S.
    // const scalar a=-0.0828;
    // const scalar b=173.59;
    // const scalar c=50.;
    // return max(a*T + b, c);
}


Foam::scalar Foam::mazarsDamageModel::epsilonEqt0
(
    const scalar& T
)
{
    return 0.00001;
}


Foam::scalar Foam::mazarsDamageModel::epsilonEqc0
(
    const scalar& T
)
{
    // Parameters by M. Reymond
    // if(T < 1373)
    // {
        return 2e8/200e9;
    // }
    // else
    // {
    //     //NO MODIF
    //     // scalar epsc = (- 191582.012847966*T + 465864989.293362) / 200e9;
    //     // epsc = (epsc < 25e6/200e9) ?  25e6/200e9 : epsc;

    //     //MODIF 1
    //     // scalar epsc = (63.4991009344333*pow(T,2) - 390707.837861043*T + 618111540.809438)/200e9;
    //     // return epsc; 

    //     //MODIF 2
    //     scalar epsc = (60.4451853901591*pow(T,2) - 362183.264348537*T + 581015887.154)/200e9;
    //     epsc = (epsc < 40e6/200e9) ?  40e6/200e9 : epsc;
    //     return epsc; 

    // }
    // //MODIF 1
    // // if(T>=2773)
    // // {
    // //     return (63.4991009344333*pow(2773,2) - 390707.837861043*2773 + 618111540.809438)/200e9;
    // // }

    // // // Parameters by JS
    // // const scalar a = 1.e6;
    // // const scalar e = -67.3409;
    // // const scalar f = 3.72591e5;
    // // const scalar g = -5.76021e8;
    // // if (T < 1500)
    // // {
    // //     return a/200e9;
    // // }
    // // if(T >= 1500 && T < 22773) 
    // // {
    // //     return mag((e*T*T+f*T+g)/200e9);
    // // }
    // // if(T >= 2773)
    // // {
    // //     return 60e6/200e9;
    // // }
}


Foam::scalar Foam::mazarsDamageModel::triaxialFactor
(
    const scalar& s1,
    const scalar& s2,
    const scalar& s3
)
{
    // positive part
    const scalar pps1 = max(0, s1);
    const scalar pps2 = max(0, s2);
    const scalar pps3 = max(0, s3);

    const scalar sn = max(mag(s1), max(mag(s2), mag(s3)));
    const scalar seps = 1e-12 * 190e9; // check stress state of the material
                                    // (190e9 --> 190 gigapascal, UO2 young modulus)
    const scalar r = (sn > seps) ? 
                    min((pps1+pps2+pps3) / (mag(s1)+mag(s2)+mag(s3)), 1)
                    : 0 ;
    return r;
}


Foam::scalar Foam::mazarsDamageModel::A
(
    const scalar& r,
    const scalar& At,
    const scalar& Ac,
    const scalar& k
)
{
    const scalar r2 = r*r;
    return At * (2*r2 * (1 - 2*k) - r*(1 - 4*k)) +  
           Ac * (2*r2 - 3*r + 1);
}


Foam::scalar Foam::mazarsDamageModel::B
(
    const scalar& r,
    const scalar& Bt,
    const scalar& Bc
)
{
    const scalar r2 = r*r;
    return Bt * pow(r, r2 - 2*r +2) +
           Bc * (1 - pow(r, r2 - 2*r + 2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::mazarsDamageModel::mazarsDamageModel
(
    const fvMesh& mesh,
    const dictionary& damageDict,
    const dictionary& baseDict
)
:
    damageModel(mesh, damageDict, baseDict),    
    d_(createOrLookup<scalar>(mesh, "damage", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    Yt_(createOrLookup<scalar>(mesh, "mazars::Yt", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    Yc_(createOrLookup<scalar>(mesh, "mazars::Yc", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    Y_(createOrLookup<scalar>(mesh, "mazars::Y", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    r_(createOrLookup<scalar>(mesh, "triaxialFactor", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    k_(1.0),
    d_max_(damageDict.lookupOrDefault("damageMax", 0.99)),
    temperatureThreshold_(
        damageDict.lookupOrDefault("temperatureFragmentationThreshold", 1500)),
    equivalentStrainRateThreshold_(
        damageDict.lookupOrDefault("equivalentStrainRateThreshold", false)),
    strainRateThreshold_(damageDict.lookupOrDefault("strainRateThreshold", 0.1)),
    excludeAxialCmpt_(damageDict.lookupOrDefault<bool>("excludeAxialCmpt", false)),
    relax_(damageDict.lookupOrDefault<scalar>("relaxDamage", 1.0)),
    currentTime_(mesh_.time().value())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mazarsDamageModel::~mazarsDamageModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::mazarsDamageModel::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev,
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    // Material properties
    const scalarField& nu(mesh_.lookupObject<volScalarField>("nu"));
    const scalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const scalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));

    // Regerence to T, epsEl and epsEl (for ddtEpsEl)
    const scalarField& T(mesh_.lookupObject<volScalarField>("T"));
    const symmTensorField& epsilonEl0(
        mesh_.lookupObject<volSymmTensorField>("epsEl").oldTime());

    // Old values of Yt and Yc
    const scalarField& Yt0(Yt_.oldTime().internalField());
    const scalarField& Yc0(Yc_.oldTime().internalField());

    if(mesh_.time().value() > currentTime_)
    {            
        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];

            // Take copy to current and old epsEl
            symmTensor epsilonEli = epsilonEl[cellI];
            symmTensor epsilonEl0i = epsilonEl0[cellI];

            // Eliminate shear strains in x-z and y-z
            if(excludeAxialCmpt_)
            {
                // Store zz cmpt
                scalar zCmpt = epsilonEli.zz();
                scalar zCmptOld = epsilonEl0i.zz();

                // Kill all -z components
                tensor R = Foam::tensor(1,0,0,0,1,0,0,0,1);
                epsilonEli = symm((R&epsilonEl[cellI])&R.T());
                epsilonEl0i = symm((R&epsilonEl0[cellI])&R.T());

                // Re-introduce zz component
                epsilonEli.zz() = zCmpt;
                epsilonEl0i.zz() = zCmptOld;
            }

            // Principal strains
            vector eigenVal = Foam::vector::zero;
            tensor eigenVec = tensor::zero;
            eig3().eigen_decomposition(epsilonEli, eigenVec, eigenVal);

            scalar& e1 = eigenVal[0];
            scalar& e2 = eigenVal[1];
            scalar& e3 = eigenVal[2];
            scalar tr = e1 + e2 + e3;

            // Principal stresses
            vector sigmaEigenVal = vector::zero;
            scalar& s1 = sigmaEigenVal.x();
            scalar& s2 = sigmaEigenVal.y();
            scalar& s3 = sigmaEigenVal.z();

            s1 = 2*mu[cellI]*e1 + lambda[cellI]*tr;
            s2 = 2*mu[cellI]*e2 + lambda[cellI]*tr;
            s3 = 2*mu[cellI]*e3 + lambda[cellI]*tr;

            // Remove axial component
            if (excludeAxialCmpt_)
            {
                int zAlignedIndex = 0;
                eigenVec = eigenVec.T();
                scalar maxZComponent = std::abs(eigenVec.vectorComponent(0)[2]);

                for (int i = 1; i < 3; ++i) 
                {
                    scalar zComponent = std::abs(eigenVec.vectorComponent(i).z());
                    if (zComponent > maxZComponent) 
                    {
                        maxZComponent = zComponent;
                        zAlignedIndex = i;
                    }
                }

                // The index `zAlignedIndex` now points to the 
                // eigenvector most aligned with the z-direction
                eigenVec = eigenVec.T();
                sigmaEigenVal[zAlignedIndex] *= 0;
                eigenVal[zAlignedIndex] *= 0;
            }

            // Triaxial factor
            r_[cellI] = triaxialFactor(s1, s2, s3);

            // mu-model Parameters
            const scalar At = this->At(T[cellI]);
            const scalar Ac = this->Ac(T[cellI]);
            const scalar Bt = this->Bt(T[cellI]);
            const scalar Bc = this->Bc(T[cellI]);
            const scalar epsilonEqt0 = this->epsilonEqt0(T[cellI]);
            const scalar epsilonEqc0 = this->epsilonEqc0(T[cellI]);

            // Calculate strain invariants
            scalar inv1 = e1 + e2 + e3;
            scalar inv2 = 0.5 * (pow(e1-e2,2) + pow(e2-e3, 2) + pow(e3-e1, 2));
            
            // Current equivalent strain - traction
            scalar epsilonEqt(inv1 > 0 ? 
                inv1/(2*(1-2*nu[cellI])) + sqrt(inv2)/(2*(1+nu[cellI]))
                : 0);

            // Update Yt
            Yt_[cellI] = max(max(Yt0[cellI], epsilonEqt0), epsilonEqt);

            // Current equivalent strain - compression
            scalar epsilonEqc(inv1 < 0 ?
                inv1/(5*(1-2*nu[cellI])) + 6*sqrt(inv2)/(5*(1+nu[cellI])) 
                : 0);

            // Update Yc
            if(equivalentStrainRateThreshold_)
            {
                // Calculate ddtEpsEl
                scalar deltaT(mesh_.time().deltaT().value());
                symmTensor ddtEpsEl((epsilonEli - epsilonEl0i)/deltaT);

                scalar ddtEpsElEq = sqrt(
                    pow(ddtEpsEl.component(symmTensor::XX), 2) 
                    + pow(ddtEpsEl.component(symmTensor::YY), 2)
                    + pow(ddtEpsEl.component(symmTensor::ZZ), 2));
                
                if(T[cellI]>temperatureThreshold_ && ddtEpsElEq >= strainRateThreshold_)
                {
                    Yc_[cellI] = max(max(Yc0[cellI], epsilonEqc0), epsilonEqc);
                }
                else
                {
                    Yc_[cellI] = max(Yc0[cellI], epsilonEqc0);
                }
            }
            else
            {
                if(T[cellI]>temperatureThreshold_)
                {
                    Yc_[cellI] = max(max(Yc0[cellI], epsilonEqc0), epsilonEqc);
                }
                else
                {
                    Yc_[cellI] = max(Yc0[cellI], epsilonEqc0);
                }
            }

            // thermodynamics variables computation
            const scalar A = this->A(r_[cellI], At, Ac, k_);
            const scalar B = this->B(r_[cellI], Bt, Bc);  
            const scalar Y0 = r_[cellI]*epsilonEqt0 + (1-r_[cellI])*epsilonEqc0;
            Y_[cellI] = r_[cellI]*Yt_[cellI] + (1-r_[cellI])*Yc_[cellI];

            // Updated damage value
            scalar d_new(
                1 - (1-A)*Y0/max(Y_[cellI], SMALL) - A*exp(-B*(Y_[cellI] - Y0))); 
            
            // Relax damage and modify E
            d_[cellI] = relax_*min(max(0, d_new), d_max_) + (1 - relax_)*d_[cellI];

            if(d_[cellI] < SMALL)
            {
                d_[cellI] = 0;
            }

            const cell& c = mesh_.cells()[cellI];

            forAll(c, faceI)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

                if (patchID > -1 and sigmaDev.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                    // References to patch fields
                    const scalarField& nup(mesh_.lookupObject<volScalarField>("nu").boundaryField()[patchID]);
                    const scalarField& mup(mesh_.lookupObject<volScalarField>("mu").boundaryField()[patchID]);
                    const scalarField& lambdap(mesh_.lookupObject<volScalarField>("lambda").boundaryField()[patchID]);

                    // Regerence to T, epsEl and epsEl (for ddtEpsEl)
                    const scalarField& Tp(mesh_.lookupObject<volScalarField>("T").boundaryField()[patchID]);
                    const symmTensorField& epsilonElp(epsilonEl.boundaryField()[patchID]);
                    const symmTensorField& epsilonEl0p(
                        mesh_.lookupObject<volSymmTensorField>("epsEl").oldTime().boundaryField()[patchID]);
                    
                    // Old values of Yt and Yc
                    const scalarField& Yt0p(Yt_.oldTime().boundaryField()[patchID]);
                    const scalarField& Yc0p(Yc_.oldTime().boundaryField()[patchID]);    

                    scalarField& rp = r_.boundaryFieldRef()[patchID];
                    scalarField& Ycp = Yc_.boundaryFieldRef()[patchID];
                    scalarField& Ytp = Yt_.boundaryFieldRef()[patchID];
                    scalarField& Yp = Y_.boundaryFieldRef()[patchID];
                    scalarField& dp = d_.boundaryFieldRef()[patchID];
                    
                    // Take copy to current and old epsEl
                    symmTensor epsilonElpf = epsilonElp[faceID];
                    symmTensor epsilonEl0pf = epsilonEl0p[faceID];   

                    // Eliminate shear strains in x-z and y-z
                    if(excludeAxialCmpt_)
                    {
                        // Store zz cmpt
                        scalar zCmpt = epsilonElpf.zz();
                        scalar zCmptOld = epsilonEl0pf.zz();

                        // Kill all -z components
                        tensor R = Foam::tensor(1,0,0,0,1,0,0,0,1);
                        epsilonElpf = symm((R&epsilonElp[faceID])&R.T());
                        epsilonEl0pf = symm((R&epsilonEl0p[faceID])&R.T());

                        // Re-introduce zz component
                        epsilonElpf.zz() = zCmpt;
                        epsilonEl0pf.zz() = zCmptOld;
                    }

                    // Principal strains
                    vector eigenVal = Foam::vector::zero;
                    tensor eigenVec = tensor::zero;
                    eig3().eigen_decomposition(epsilonElpf, eigenVec, eigenVal);

                    scalar& e1 = eigenVal[0];
                    scalar& e2 = eigenVal[1];
                    scalar& e3 = eigenVal[2];
                    scalar tr = e1 + e2 + e3;

                    // Principal stresses
                    vector sigmaEigenVal = vector::zero;
                    scalar& s1 = sigmaEigenVal.x();
                    scalar& s2 = sigmaEigenVal.y();
                    scalar& s3 = sigmaEigenVal.z();

                    s1 = 2*mup[faceID]*e1 + lambdap[faceID]*tr;
                    s2 = 2*mup[faceID]*e2 + lambdap[faceID]*tr;
                    s3 = 2*mup[faceID]*e3 + lambdap[faceID]*tr;

                    // Remove axial component
                    if (excludeAxialCmpt_)
                    {
                        int zAlignedIndex = 0;
                        eigenVec = eigenVec.T();
                        scalar maxZComponent = std::abs(eigenVec.vectorComponent(0)[2]);

                        for (int i = 1; i < 3; ++i) 
                        {
                            scalar zComponent = std::abs(eigenVec.vectorComponent(i).z());
                            if (zComponent > maxZComponent) 
                            {
                                maxZComponent = zComponent;
                                zAlignedIndex = i;
                            }
                        }

                        // The index `zAlignedIndex` now points to the 
                        // eigenvector most aligned with the z-direction
                        eigenVec = eigenVec.T();
                        sigmaEigenVal[zAlignedIndex] *= 0;
                        eigenVal[zAlignedIndex] *= 0;
                    }

                    // Triaxial factor
                    rp[faceID] = triaxialFactor(s1, s2, s3);

                    // mu-model Parameters
                    const scalar At = this->At(Tp[faceID]);
                    const scalar Ac = this->Ac(Tp[faceID]);
                    const scalar Bt = this->Bt(Tp[faceID]);
                    const scalar Bc = this->Bc(Tp[faceID]);
                    const scalar epsilonEqt0 = this->epsilonEqt0(Tp[faceID]);
                    const scalar epsilonEqc0 = this->epsilonEqc0(Tp[faceID]);
                    
                    // Calculate strain invariants
                    scalar inv1 = e1 + e2 + e3;
                    scalar inv2 = 0.5 * (pow(e1-e2,2) + pow(e2-e3, 2) + pow(e3-e1, 2));
                    
                    // Current equivalent strain - traction
                    scalar epsilonEqt(inv1 > 0 ? 
                        inv1/(2*(1-2*nup[faceID])) + sqrt(inv2)/(2*(1+nup[faceID]))
                        : 0);

                    // Update Yt
                    Ytp[faceID] = max(max(Yt0p[faceID], epsilonEqt0), epsilonEqt);

                    // Current equivalent strain - compression
                    scalar epsilonEqc(inv1 < 0 ?
                        inv1/(5*(1-2*nup[faceID])) + 6*sqrt(inv2)/(5*(1+nup[faceID])) 
                        : 0);

                    // Update Yc
                    if(equivalentStrainRateThreshold_)
                    {
                        // Calculate ddtEpsEl
                        scalar deltaT(mesh_.time().deltaT().value());
                        symmTensor ddtEpsEl((epsilonElpf - epsilonEl0pf)/deltaT);

                        scalar ddtEpsElEq = sqrt(
                            pow(ddtEpsEl.component(symmTensor::XX), 2) 
                            + pow(ddtEpsEl.component(symmTensor::YY), 2)
                            + pow(ddtEpsEl.component(symmTensor::ZZ), 2));
                        
                        if(Tp[faceID]>temperatureThreshold_ && ddtEpsElEq >= strainRateThreshold_)
                        {
                            Ycp[faceID] = max(max(Yc0p[faceID], epsilonEqc0), epsilonEqc);
                        }
                        else
                        {
                            Ycp[faceID] = max(Yc0p[faceID], epsilonEqc0);
                        }
                    }
                    else
                    {
                        if(Tp[faceID]>temperatureThreshold_)
                        {
                            Ycp[faceID] = max(max(Yc0p[faceID], epsilonEqc0), epsilonEqc);
                        }
                        else
                        {
                            Ycp[faceID] = max(Yc0p[faceID], epsilonEqc0);
                        }
                    }

                    // thermodynamics variables computation
                    const scalar A = this->A(rp[faceID], At, Ac, k_);
                    const scalar B = this->B(rp[faceID], Bt, Bc);  
                    const scalar Y0 = rp[faceID]*epsilonEqt0 + (1-rp[faceID])*epsilonEqc0;
                    Yp[faceID] = rp[faceID]*Ytp[faceID];

                    // Updated damage value
                    scalar d_new(
                        1 - (1-A)*Y0/max(Yp[faceID], SMALL) - A*exp(-B*(Yp[faceID] - Y0))); 
                    
                    // Relax damage and modify E
                    dp[faceID] = relax_*min(max(0, d_new), d_max_) + (1 - relax_)*dp[faceID];

                    if(dp[faceID] < SMALL)
                    {
                        dp[faceID] = 0;
                    }
                }
            }
        }
    }
 
    // Applying damage to sigmaDev and sigmaHyd
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        if (excludeAxialCmpt_)
        {
            volTensorField& dTensor = 
            (createOrLookup<tensor>(mesh_, "damageTensor", dimless, Foam::tensor::zero, 
                calculatedFvPatchField<tensor>::typeName));
            
            tensor& D = dTensor.ref()[cellI];
            D = d_[cellI]*Foam::tensor(1,1,1,1,1,1,1,1,0);

            // Build sigma and apply damage
            symmTensor sigma = (sigmaHyd[cellI]*I) + sigmaDev[cellI];
            sigma = Foam::cmptMultiply(symm(Foam::tensor::one - D), sigma);

            // Rebuild sigmaHyd and sigmaDev
            sigmaHyd[cellI] = tr(sigma)/3.0;
            sigmaDev[cellI] = dev(sigma);
        }
        else
        {
            sigmaHyd[cellI] *= (1 - d_[cellI]);
            sigmaDev[cellI] *= (1 - d_[cellI]);
        }

        const cell& c = mesh_.cells()[cellI];

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if (patchID > -1 and sigmaDev.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                
                scalarField& sigmaHydF = 
                sigmaHyd.boundaryFieldRef()[patchID];

                symmTensorField& sigmaDevF = 
                sigmaDev.boundaryFieldRef()[patchID];

                const scalarField& dF = 
                d_.boundaryField()[patchID];

                if (excludeAxialCmpt_)
                {
                    volTensorField& dTensor = 
                    (createOrLookup<tensor>(mesh_, "damageTensor", dimless, Foam::tensor::zero, 
                        calculatedFvPatchField<tensor>::typeName));
            
                    tensor& D = dTensor.boundaryFieldRef()[patchID][faceID];
                    D = dF[faceID]*Foam::tensor(1,1,1,1,1,1,1,1,0);

                    // Build sigma and apply damage
                    symmTensor sigma = (sigmaHydF[faceID]*I) + sigmaDevF[faceID];
                    sigma = Foam::cmptMultiply(symm(Foam::tensor::one - D), sigma);

                    // Rebuild sigmaHyd and sigmaDev
                    sigmaHydF[faceID] = tr(sigma)/3.0;
                    sigmaDevF[faceID] = dev(sigma);
                }
                else
                {
                    sigmaHydF[faceID] *= (1 - dF[faceID]);
                    sigmaDevF[faceID] *= (1 - dF[faceID]);
                }
            }
        }
    }
}

// ************************************************************************* 