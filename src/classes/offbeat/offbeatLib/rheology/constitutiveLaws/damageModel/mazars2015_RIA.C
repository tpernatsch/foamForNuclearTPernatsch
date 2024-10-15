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

#include "mazars2015_RIA.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
#include "Time.H"
#include "eig3.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mazars2015_RIA, 0);
    addToRunTimeSelectionTable
    (
        damageModel, 
        mazars2015_RIA, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::mazars2015_RIA::At
(
    const scalar& T
)
{
    return 3.5;
}


Foam::scalar Foam::mazars2015_RIA::Ac
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


Foam::scalar Foam::mazars2015_RIA::Bt
(
    const scalar& T
)
{
    return 1750;
}


Foam::scalar Foam::mazars2015_RIA::Bc
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


Foam::scalar Foam::mazars2015_RIA::epsilonEqt0
(
    const scalar& T
)
{
    return 0.00001;
}


Foam::scalar Foam::mazars2015_RIA::epsilonEqc0
(
    const scalar& T
)
{
    // Parameters by M. Reymond
    if(T < 1373)
    {
        return 2e8/200e9;
    }
    else
    {
        //NO MODIF
        // scalar epsc = (- 191582.012847966*T + 465864989.293362) / 200e9;
        // epsc = (epsc < 25e6/200e9) ?  25e6/200e9 : epsc;

        //MODIF 1
        // scalar epsc = (63.4991009344333*pow(T,2) - 390707.837861043*T + 618111540.809438)/200e9;
        // return epsc; 

        //MODIF 2
        scalar epsc = (60.4451853901591*pow(T,2) - 362183.264348537*T + 581015887.154)/200e9;
        epsc = (epsc < 40e6/200e9) ?  40e6/200e9 : epsc;
        return epsc; 

    }
    //MODIF 1
    // if(T>=2773)
    // {
    //     return (63.4991009344333*pow(2773,2) - 390707.837861043*2773 + 618111540.809438)/200e9;
    // }

    // // Parameters by JS
    // const scalar a = 1.e6;
    // const scalar e = -67.3409;
    // const scalar f = 3.72591e5;
    // const scalar g = -5.76021e8;
    // if (T < 1500)
    // {
    //     return a/200e9;
    // }
    // if(T >= 1500 && T < 22773) 
    // {
    //     return mag((e*T*T+f*T+g)/200e9);
    // }
    // if(T >= 2773)
    // {
    //     return 60e6/200e9;
    // }
}


Foam::scalar Foam::mazars2015_RIA::triaxialFactor
(
    const symmTensor& epsilon,
    const scalar& mu,
    const scalar& lambda,
    const scalar& tr
)
{
    tensor eigenVectors = tensor::zero;
    vector eigenVal = vector::zero;

    // Compute eigen values of the strain tensor
    eig3().eigen_decomposition(epsilon, eigenVectors, eigenVal);
    scalar& e1 = eigenVal.component(0);
    scalar& e2 = eigenVal.component(1);
    scalar& e3 = eigenVal.component(2);

    // Compute principal stresses 
    const scalar s1 = 2*mu*e1 + lambda*tr;
    const scalar s2 = 2*mu*e2 + lambda*tr;
    const scalar s3 = 2*mu*e3 + lambda*tr;

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


Foam::scalar Foam::mazars2015_RIA::A
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


Foam::scalar Foam::mazars2015_RIA::B
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
Foam::mazars2015_RIA::mazars2015_RIA
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    damageModel(mesh, lawDict),    
    d_(createOrLookup<scalar>(mesh, "damage", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    Yt_(createOrLookup<scalar>(mesh, "Yt", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    Yc_(createOrLookup<scalar>(mesh, "Yc", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    Y_(createOrLookup<scalar>(mesh, "Y", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    r_(createOrLookup<scalar>(mesh, "triaxialFactor", dimless, 0.0, 
        zeroGradientFvPatchField<scalar>::typeName)),
    k_(1.0),
    d_max_(lawDict.lookupOrDefault("dmax", 0.99)),
    temperatureThreshold_(
        lawDict.lookupOrDefault("temperatureFragmentationThreshold", 1500)),
    equivalentStrainRateThreshold_(
        lawDict.lookupOrDefault("equivalentStrainRateThreshold", false)),
    strainRateThreshold_(lawDict.lookupOrDefault("strainRateThreshold", 0.1)),
    relax_(lawDict.lookupOrDefault<scalar>("relaxDamage", 1.0)),
    currentTime_(mesh_.time().value())   
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mazars2015_RIA::~mazars2015_RIA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::mazars2015_RIA::correct
(
    const labelList& addr
)
{
    // Material properties
    const scalarField& nu(mesh_.lookupObject<volScalarField>("nu"));
    const scalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const scalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));

    // Regerence to T, epsEl and epsEl (for ddtEpsEl)
    const scalarField& T(mesh_.lookupObject<volScalarField>("T"));    
    const symmTensorField& epsilonEl(
        mesh_.lookupObject<volSymmTensorField>("epsEl"));
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

            // mu-model Parameters
            const scalar At = this->At(T[cellI]);
            const scalar Ac = this->Ac(T[cellI]);
            const scalar Bt = this->Bt(T[cellI]);
            const scalar Bc = this->Bc(T[cellI]);
            const scalar epsilonEqt0 = this->epsilonEqt0(T[cellI]);
            const scalar epsilonEqc0 = this->epsilonEqc0(T[cellI]);

            // Triaxial factor calculation
            r_[cellI] = triaxialFactor(epsilonEl[cellI], mu[cellI], 
                                            lambda[cellI], tr(epsilonEl[cellI]));
            
            // Calculate strain eigenvalues
            tensor eigenVectors = tensor::zero;
            vector eigenVal = vector::zero;

            eig3().eigen_decomposition(epsilonEl[cellI], eigenVectors, eigenVal);
            scalar& e1(eigenVal.component(0));
            scalar& e2(eigenVal.component(1));
            scalar& e3(eigenVal.component(2));

            // 1st & 2nd invariant 
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
                symmTensor ddtEpsEl((epsilonEl[cellI] - epsilonEl0[cellI])/deltaT);

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
        }
    }
}

// ************************************************************************* 