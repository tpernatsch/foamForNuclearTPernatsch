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

#include "poreVelocityMOXVaporPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreVelocityMOXVaporPressure, 0);
    addToRunTimeSelectionTable(poreVelocityModel, poreVelocityMOXVaporPressure, dictionary);

}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::scalar Foam::poreVelocityMOXVaporPressure::calcPartialPressureO2
(
    const scalar OMi, 
    const scalar Ti
)
{
    //  Useful quantities
    const scalar x  = 2-OMi;
    const scalar R  = 8.3145;
    const scalar RT = R*Ti;

    // Compute equilibrium constants
    const scalar cPu    = wPu_ + 2.5*wAm_;
    const scalar kvo_ki = exp((44.0+55.8*cPu)/R)*exp(-3.76e5/RT);
    const scalar koi_ki = exp((-22.8-84.5*cPu)/R)*exp(1.05e5/RT);
    const scalar kre    = exp((153.5-96.5*cPu+331.0*pow(cPu,2.0))/R)*exp(-8.91e5/RT);

    scalar pO2(0);

    if ( OMi == 2 )
    {
        pO2 = kvo_ki/koi_ki;
    }
    else if ( OMi > 2)
    {
        // TODO: For hyperstoichiometric MOX fuels this correlation is wrong. 
        //       A better approach would be numerically solving the equation of 
        //       OM ratio given in DOI: 10.1016/j.jnucmat.2017.01.056
        pO2 = pow(koi_ki*(2-OMi), 2.0)+kvo_ki/koi_ki;
    }
    else
    {
        scalar kren4 = exp((68.8+131.3*cPu)/R)*exp(-5.15e5/RT);
        scalar kren3 = kre;

        // a,b,c term are capped in order to produce a max result of 1e300. This
        // is done because of problems arising at low temperature, where pO2 
        // tend to 0 atm. The tiny number computed by the correlation for pO2 
        // causes overflow.
        scalar a = pow(min((2-OMi)/(kvo_ki), 1e30), 10.0);
        scalar b = pow(min((2-OMi)/(pow(kren4,0.5)),1e15) ,20.0);
        scalar c = pow(min((2-OMi)/(pow(2*kren3,1.0/3.0)), 1e20),15.0);
        scalar d = 1.0/(((2-x)-(2-cPu/2.0))*80);
        
        scalar F = Foam::log10(a + b + c) + d;

        pO2 = pow(pow(10,F/5.0)+koi_ki/kvo_ki, -1.0);
    }

    return pO2;
}

const Foam::scalar Foam::poreVelocityMOXVaporPressure::calcIntegral
(
    const scalar upperLimit,
    const scalar Ti,
    const scalar q
)
{
    // Define number of integration steps
    // const scalar n = 100;
    const scalar n = 10;
    const scalar dm1 = upperLimit / n;

    scalar sum = 0;

    for (scalar i = 1; i < n; i++)
    {
        scalar m1_avg = 0.5 * (i * dm1 + (i - 1) * dm1);
        scalar devFromStoch = m1_avg * q; 
        scalar oxyMetRatio = 2 - devFromStoch;
        sum += log(calcPartialPressureO2(oxyMetRatio, Ti)) * dm1; 
    }

    return sum;
}

Foam::scalar Foam::poreVelocityMOXVaporPressure::calcVaporPressure
(
    const scalar OMi,
    const scalar Ti
)
{
    // Define useful quantities
    scalar q(wPu_);
    scalar r(wAm_);
    scalar x(2-OMi);
    scalar m(x/q);
    scalar n = (wAm_ > 0.0) ? (x/r) : 0;
    scalar R(8.3145);
    scalar RT(R*Ti);

    // Compute partial pressure O2
    const scalar pO2 = calcPartialPressureO2(OMi, Ti);

    // Gibbs energies
    scalar DG_UO2vap    = 5.67e5 - 150 * Ti;
    scalar DG_PuO2vap   = 5.71e5 - 150 * Ti;
    scalar DG_AmO2vap   = -RT*log(pow(10, 7.28-28260/Ti));

    scalar DG_A         = -0.5*RT*calcIntegral(m, Ti, q) + DG_PuO2vap;
    scalar DG_B         = (wAm_ > 0.0)
    ? -0.5*RT*calcIntegral(n, Ti, r) + DG_AmO2vap
    : 0;
    
    scalar DG_UOUO2     = -4.71e5 + 71 * Ti;
    scalar DG_PuOPuO2   = -3.52e5 + 69 * Ti;
    
    scalar DG_UUO       = -5.28e5 + 62 * Ti;
    scalar DG_PuPuO     = -4.98e5 + 46 * Ti;
    
    scalar DG_UO2UO3    = -4.04e5 + 90 * Ti;
    
    scalar DG_AmO2Am    = (-3.883e5 + 30.7 * Ti ) - (263650 - 112.32 * Ti );
    scalar DG_AmOAm     = (-1.377e5 - 28.8 * Ti ) - (263650 - 112.32 * Ti );

    // Compute partial pressures

    scalar p_UO2    = (1-q-r)*exp(-DG_UO2vap/RT);
    scalar p_PuO2   = (q*pow(pO2,m/2.0))*exp(-DG_A/RT);
    scalar p_AmO2   = r*pow(pO2,n/2.0)*exp(-DG_B/RT);
    
    scalar p_UO     = p_UO2/(exp(-DG_UOUO2/RT)*pow(pO2,0.5));
    scalar p_PuO    = p_PuO2/(exp(-DG_PuOPuO2/RT)*pow(pO2,0.5));

    scalar p_U      = p_UO/(exp(-DG_UUO/RT)*pow(pO2,0.5));
    scalar p_Pu     = p_PuO/(exp(-DG_PuPuO/RT)*pow(pO2,0.5));
    
    scalar p_UO3    = exp(-DG_UO2UO3/RT)*p_UO2*pow(pO2,0.5);
    
    // scalar p_AmO    = p_AmO2/(exp(-DG_AmO2Am/RT)*pow(pO2,0.5));
    // scalar p_Am     = p_AmO/(exp(-DG_AmOAm/RT)*pow(pO2,0.5));

    scalar p_AmO    = p_AmO2/(exp(-DG_AmOAm/RT)*pow(pO2,0.5));
    scalar p_Am     = p_AmO/(exp(-DG_AmO2Am/RT)*pow(pO2,0.5));

    //  Compute total vapor pressure in Pascal
    scalar atm2pa(101325);

    return atm2pa * (p_UO2+p_PuO2+p_UO+p_U+p_UO3+p_PuO+p_Pu+p_AmO2+p_Am+p_AmO);

}

Foam::scalar Foam::poreVelocityMOXVaporPressure::calcVaporPressureDerivative
(
    const scalar OMi,
    const scalar Ti
)
{
    const scalar P1 = calcVaporPressure(OMi, Ti-1);
    const scalar P2 = calcVaporPressure(OMi, Ti+1);

    return (P2-P1)/2.0;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreVelocityMOXVaporPressure::poreVelocityMOXVaporPressure
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
:
    poreVelocityModel(mesh, dict, defaultModel),
    T_(nullptr),
    gradT_(nullptr),
    OM_(nullptr),
    CPuO2_(5.396e-8),
    CUO2_(5.4704e-8),
    wPu_(dict.lookup<scalar>("ratioPuMetal")),
    wU_(dict.lookup<scalar>("ratioUMetal")),
    wAm_(dict.lookup<scalar>("ratioAmMetal")),
    xPuO2_((wPu_ / 276)/(wU_ / 270 + wPu_ / 276 + wAm_ / 275)),
    xUO2_((wU_ / 270)/(wU_ / 270 + wPu_ / 276  + wAm_ / 275)),
    kb_(1.3807e-16),
    d_(4.33e-8),
    m_(28+270),
    nD_(1.4133548e6)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poreVelocityMOXVaporPressure::~poreVelocityMOXVaporPressure()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::poreVelocityMOXVaporPressure::correct
(
    const labelList& addr
)
{
    if( T_ == nullptr )
    {
        T_ = &mesh_.lookupObject<volScalarField>("T");
    }

    if( gradT_ == nullptr )
    {
        gradT_ = &mesh_.lookupObject<volVectorField>("gradT");
    }

    if( OM_ == nullptr )
    {
        OM_ = &mesh_.lookupObject<volScalarField>("oxygenMetalRatio");
    }

    // Reference to pore velocity field
    vectorField& poreVel = poreVelocity_.ref();

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T_->internalField()[cellI];

        const scalar OMi = OM_->internalField()[cellI];

        // PI constant
        const scalar PI = 3.1459;
        
        const vector gradTi = gradT_->internalField()[cellI];

        // Compute fuel molecular volume in cm3
        scalar omega = 0.25*pow(xPuO2_*CPuO2_ + xUO2_*CUO2_, 3.0);

        // Compute total vapor pressure
        scalar P = calcVaporPressure(OMi, Ti);

        // Compute dP/dT
        scalar dPdT = calcVaporPressureDerivative(OMi, Ti);

        // Convert P and dP/dT to suitable unit of measure (g, cm)
        P    *= 10;
        dPdT *= 10;

        // Compute number of gas molecules per unit vol (ideal gas law)
        scalar n = P/(kb_*Ti);

        // Compute diffusion coeff of PuO2 and UO2 gas in the void
        scalar D12 = 3.0/(8.0*PI*nD_*pow(d_,2))*sqrt((PI*kb_*Ti)/(2*m_));

        // Compute dn/dT
        const scalar dndT = 1.0/(kb_*pow(Ti,2.0))*(dPdT*Ti-P);

        // Pore velocity in (cm/s)
        poreVel[cellI] =  omega * D12 * dndT * gradTi/100.0 * alphaPoreVelocity_;

        // Convert pore velocity in m/s
        poreVel[cellI] /= 100.0;

        // 
        // Correct boundary values
        // 

        // Take reference of current cell
        const cell& c = mesh_.cells()[cellI];  

        // Loop over all faces of current cell
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            const polyPatch& pp(mesh_.boundaryMesh()[patchID]);

                    
            if (patchID > -1 and poreVelocity_.boundaryField()[patchID].size())
            {
                const vectorField deltas = pp.faceCentres() - pp.faceCellCentres();
                
                const label faceID = pp.whichFace(c[faceI]);

                // Take references boundary T
                const scalarField& Tp(T_->boundaryField()[patchID]);

                // Compute grad T at the boundary
                const vector delta = deltas[faceID];

                const scalar magGradTp = (Tp[faceID]-Ti)/mag(delta);
                const vector gradTp = magGradTp * delta/mag(delta);

                // Take references boundary OM
                const scalarField& OMp(OM_->boundaryField()[patchID]);
                
                // Take references boundary v
                vectorField& vP = poreVelocity_.boundaryFieldRef()[patchID];

                // Compute total vapor pressure
                scalar Pp = calcVaporPressure(OMp[faceID], Tp[faceID]);

                // Compute dP/dT
                scalar dPdTp = calcVaporPressureDerivative(OMp[faceID], Tp[faceID]);

                // Convert P and dP/dT to suitable unit of measure (g, cm)
                Pp    *= 10;
                dPdTp *= 10;

                // Compute number of gas molecules per unit vol (ideal gas law)
                scalar np = Pp/(kb_*Tp[faceID]);

                // Compute diffusion coeff of PuO2 and UO2 gas in the void
                scalar D12p = 3.0/(8.0*PI*nD_*pow(d_,2))*sqrt((PI*kb_*Tp[faceID])/(2*m_));

                // Compute dn/dT
                const scalar dndTp = 1.0/(kb_*pow(Tp[faceID],2.0))*(dPdTp*Tp[faceID]-Pp);
                
                // Compute velocity in cm/s
                vP[faceID] = omega * D12p * dndTp * gradTp/100 * alphaPoreVelocity_;

                // Convert velocity to m/s
                vP[faceID] /= 100.0;

            }
        }
    }
}

// ************************************************************************* //
