/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "nuclearFuelPin.H"
#include "addToRunTimeSelectionTable.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(nuclearFuelPin, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        nuclearFuelPin, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelPin::nuclearFuelPin
(
    const structureModel& structure,
    const dictionary& dict,
    const word name
)
:
    powerModel
    (
        structure,
        dict,
        name
    ),
    Trad_
    (
        IOobject
        (
            IOobject::groupName(typeName, "Trad."+name),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_.cells().size()
    ),
    powerDensity_
    (
        mesh_.lookupObjectRef<volScalarField>("nuclearFuelPowerDensity")
    ),
    TavFuel_
    (
        mesh_.lookupObjectRef<volScalarField>("TavFuel")
    ),
    TavClad_
    (
        mesh_.lookupObjectRef<volScalarField>("TavClad")
    ),
    Tfi_
    (
        IOobject
        (
            "Tfi."+typeName+"."+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfo_
    (
        IOobject
        (
            "Tfo."+typeName+"."+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tci_
    (
        IOobject
        (
            "Tci."+typeName+"."+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tco_
    (
        IOobject
        (
            "Tco."+typeName+"."+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    fuelSubMeshSize_(this->lookupType<label>("fuelMeshSize")),
    cladSubMeshSize_(this->lookupType<label>("cladMeshSize")),
    subMeshSize_(fuelSubMeshSize_+cladSubMeshSize_),
    r_(0),
    rfi_(this->lookupType<scalar>("fuelInnerRadius")),
    rfo_(this->lookupType<scalar>("fuelOuterRadius")),
    rci_(this->lookupType<scalar>("cladInnerRadius")),
    rco_(this->lookupType<scalar>("cladOuterRadius")),
    drf_((rfo_-rfi_)/(fuelSubMeshSize_-1)),
    drc_((rco_-rci_)/(cladSubMeshSize_-1)),
    drg_(rci_-rfo_),
    rhof_(this->lookupType<scalar>("fuelDensity")),
    rhoc_(this->lookupType<scalar>("cladDensity")),
    kf_(this->lookupType<scalar>("fuelConductivity")),
    kc_(this->lookupType<scalar>("cladConductivity")),
    Cpf_(this->lookupType<scalar>("fuelHeatCapacity")),
    Cpc_(this->lookupType<scalar>("cladHeatCapacity")),
    gapH_(this->lookupType<scalar>("gapConductance")),
    hollowFuel_( (rfi_ < VSMALL) ? true : false)
{   
    //- Init radial mesh
    r_.append(rfi_);
    for (int i = 0; i < fuelSubMeshSize_-1; i++)
    {
        r_.append(r_.last() + drf_);
    }
    r_.append(r_.last() + drg_);
    for (int i = 0; i < cladSubMeshSize_-1; i++)
    {
        r_.append(r_.last() + drc_);
    }

    //- Init empty Trad
    forAll(mesh_.cells(), i)
    {
        Trad_.set(i, new Field<scalar>(0, 0));
    }

    //- If the power keyword is found, write it to the powerDensity field
    if (this->found("powerDensity"))
    {
        scalar q(this->lookupType<scalar>("powerDensity"));
        forAll(this->cellList_, i)
        {
            label celli(this->cellList_[i]);
            powerDensity_[celli] = q;
        }
        powerDensity_.correctBoundaryConditions();
    }

    //- If the files are present, reconstruct initial Trad_ profile 
    //  analytically. The analytical form is:
    //  T(r) = -(1/4)*powerDensity_(r)*r^2/k + C*ln(r)/k + D
    //  with C and D coming from imposing fixedValue BC on all sides,
    //  equal to the starting temperatures found in the files.
    if 
    (
            Tfi_.typeHeaderOk<volScalarField>(true)
        and Tfo_.typeHeaderOk<volScalarField>(true)
        and Tci_.typeHeaderOk<volScalarField>(true)
        and Tco_.typeHeaderOk<volScalarField>(true)
    )
    {
        Info<< "Found pin " << name <<" temperatures"
            <<" reconstructing profiles " << endl;
        forAll(this->cellList_, i)
        {
            label celli(this->cellList_[i]);
            Trad_.set(celli, new Field<scalar>(subMeshSize_, 0));
            scalar q = powerDensity_[celli];
            scalar tfi = Tfi_[celli];
            scalar tfo = Tfo_[celli];
            scalar tci = Tci_[celli];
            scalar tco = Tco_[celli];
            scalar Cf;
            scalar Df;
            scalar Cc;
            scalar Dc;
            if (hollowFuel_)
            {
                Cf = 0;
            }
            else
            {
                Cf = 
                    (kf_*(tfi-tfo)-0.25*q*(sqr(rfo_)-sqr(rfi_)))
                    /(log(rfi_/rfo_));
            }
            Df = tfo+(0.25*q*sqr(rfo_)-Cf*log(rfo_))/kf_;
            Cc = (tco-tci)*kc_/(log(rco_/rci_));
            Dc = tci - log(rci_)*Cc/kc_;
            
            forAll(Trad_[celli], j)
            {
                scalar r(r_[j]);
                if (j < fuelSubMeshSize_)
                {
                    Trad_[celli][j] = -0.25*q*sqr(r)/kf_ + Df;
                    Trad_[celli][j] +=
                        (hollowFuel_) ? (0.0) : (Cf*log(r)/kf_);
                }
                else
                {   
                    Trad_[celli][j] = Cc*log(r)/kc_ + Dc;
                }
            }
            
            Info << Trad_[celli] << endl;
        }
    }
    else //- Otherwise, read from dict
    {
        Info<< "Reading initial temperatures for pin " << name
            << " from dictionary"<< endl;
        scalar Tf0(this->lookupType<scalar>("fuelT"));
        scalar Tc0(this->lookupType<scalar>("cladT"));
        forAll(this->cellList_, i)
        {
            label celli(this->cellList_[i]);
            Trad_.set(celli, new Field<scalar>(subMeshSize_, 0));
            forAll(Trad_[celli], subCelli)
            {
                Trad_[celli][subCelli] = 
                    (subCelli < fuelSubMeshSize_) ?
                    Tf0 : Tc0;
            }
        }
    }

    //- Set I/O fields
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        Tfi_[celli] = Trad_[celli][0];
        Tfo_[celli] = Trad_[celli][fuelSubMeshSize_-1];
        Tci_[celli] = Trad_[celli][fuelSubMeshSize_];
        Tco_[celli] = Trad_[celli][subMeshSize_-1]; 
    }
    Tfi_.correctBoundaryConditions();
    Tfo_.correctBoundaryConditions();
    Tci_.correctBoundaryConditions();
    Tco_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::nuclearFuelPin::~nuclearFuelPin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::powerModels::nuclearFuelPin::updateLocalTemperatureProfile
(
    label celli,
    scalar HTSumi,
    scalar HSumi
)
{
    //- Init matrix, source
    SquareMatrix<scalar> M(subMeshSize_, subMeshSize_, Foam::zero());
    List<scalar> S(subMeshSize_, 0.0);

    //- Recurrent quantities
    scalar irf2(1.0/sqr(drf_));
    scalar irc2(1.0/sqr(drc_));

    scalar q(powerDensity_[celli]);

    scalar Xf(rhof_*Cpf_/mesh_.time().deltaT().value());
    scalar Xc(rhoc_*Cpc_/mesh_.time().deltaT().value());
    Field<scalar>& TOld = Trad_.oldTime()[celli];

    //- Construct matrix and source
    {
        //- Set zeroGradient BC at fuel inner surface
        M[0][1] =   -2*irf2*kf_;
        M[0][0] =   -M[0][1]+Xf;
        S[0] =      q+TOld[0]*Xf;

        //- Fuel bulk
        for (int i = 1; i < fuelSubMeshSize_-1; i++)
        {
            scalar r(r_[i]);

            M[i][i+1] =     -kf_*(irf2+1.0/(2*r*drf_));
            M[i][i-1] =     -kf_*(irf2-1.0/(2*r*drf_));
            M[i][i] =       -M[i][i+1]-M[i][i-1]+Xf;
            S[i] =          q+TOld[i]*Xf;
        }

        //- Fuel outer surface, convective BC with inner cladding surface via
        //  gap conductance
        {
            label i(fuelSubMeshSize_-1);
            scalar A(kf_*irf2);
            scalar B(kf_/(2*rfo_*drf_));
            scalar C(2*drf_/kf_);

            M[i][i+1] =     -gapH_*C*(A+B);
            M[i][i-1] =     -2*A;
            M[i][i] =       -M[i][i+1]-M[i][i-1]+Xf;
            S[i] =          q+TOld[i]*Xf;
        }
        
        //- Cladding inner surface, convective BC with outer fuel surface via
        //  gap conductance adjusted by radii ratio (to preserve total heat
        //  flow as geometry is cylindrical)
        //- EDIT: while the corrective ratio should be (rfo_/rci_), the actual
        //  factor that works is (rfo_/(rci_+drc_)). I have no mathematical 
        //  justification for this, I simply had an intuition that it just
        //  might help avoid temperature profile overestimation for coarse
        //  pin meshes (for fine meshes drc_ -> 0 and it does not make a
        //  difference) or for very low fluid heat transfer coefficient. And
        //  it does work!
        {
            label i(fuelSubMeshSize_);
            scalar A(kc_*irc2);
            scalar B(kc_/(2*rci_*drc_));
            scalar C(2*drc_*(rfo_/(rci_+drc_))/kc_);

            M[i][i-1] =     -gapH_*C*(A+B);
            M[i][i+1] =     -2*A;
            M[i][i] =       -M[i][i+1]-M[i][i-1]+Xc;
            S[i] =          TOld[i]*Xc;
        }

        //- Cladding bulk
        for (int i = fuelSubMeshSize_+1; i < subMeshSize_-1; i++)
        {
            scalar r(r_[i]);
            M[i][i+1] = -kc_*(irc2+1.0/(2*r*drc_));
            M[i][i-1] = -kc_*(irc2-1.0/(2*r*drc_));
            M[i][i] =   -M[i][i+1]-M[i][i-1]+Xc;
            S[i] =      TOld[i]*Xc;
        }

        //- Cladding outer surface, convective BC with fluid(s) wetting the pin
        {
            label i(subMeshSize_-1);
            scalar A(kc_*irc2);
            scalar B(kc_/(2*rci_*drc_));
            scalar C(2*drc_/kc_);

            M[i][i-1] =     -2*A;
            M[i][i] =       2*A + HSumi*C*(A+B)+Xc;
            S[i] =          HTSumi*C*(A+B)+TOld[i]*Xc;
        }
    }
    
    //- Solve linear system
    solve(Trad_[celli], M, S);

    //- Set fields (inner/outer fuel/clad)
    Tfi_[celli] = Trad_[celli][0];
    Tfo_[celli] = Trad_[celli][fuelSubMeshSize_-1];
    Tci_[celli] = Trad_[celli][fuelSubMeshSize_];
    Tco_[celli] = Trad_[celli][subMeshSize_-1];

    /* Check adjusted flux conservation
    scalar gradQc(kc_*(Trad_[celli][subMeshSize_-2]-Tco_[celli])/drc_);
    scalar gradQh(gapH_*(Tfo_[celli]-Tci_[celli])*(rfo_)/(rco_));
    Info << gradQc << " " << gradQh << endl;
    */

    //- Compute and set average temperatures
    
    //- Fuel
    scalar totAreaf(sqr(rfo_)-sqr(rfi_));
    scalar Tavf(0);
    for (int i = 0; i < fuelSubMeshSize_; i++)
    {
        scalar ri(r_[i]);
        Tavf += 
                Trad_[celli][i]
            *   drf_*(2*ri+drf_)/totAreaf;
    }
    TavFuel_[celli] = Tavf;
    
    //- Cladding
    scalar totAreac(sqr(rco_)-sqr(rci_));
    scalar Tavc(0);
    for (int i = fuelSubMeshSize_; i < subMeshSize_; i++)
    {
        scalar ri(r_[i]);
        Tavc += 
                Trad_[celli][i]
            *   drc_*(2*ri+drc_)/totAreac;
    }
    TavClad_[celli] = Tavc;
}


Foam::tmp<Foam::volScalarField> Foam::powerModels::nuclearFuelPin::T() const
{
    tmp<volScalarField> tT
    (
        new volScalarField
        (
            Tco_
        )
    );
    volScalarField& T = tT.ref();
    T.primitiveFieldRef() *= cellField_;
    T.correctBoundaryConditions();
    return tT;
}


void Foam::powerModels::nuclearFuelPin::correct
(
    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
)
{
    forAll(this->cellList_, i)
    {
        label celli(this->cellList_[i]);
        updateLocalTemperatureProfile(celli, HTSum[celli], HSum[celli]);
    }
    Tfi_.correctBoundaryConditions();
    Tfo_.correctBoundaryConditions();
    Tci_.correctBoundaryConditions();
    Tco_.correctBoundaryConditions();
    TavFuel_.correctBoundaryConditions();
    TavClad_.correctBoundaryConditions();
}


// ************************************************************************* //


/*-------------------------- EQUATION DERIVATION ----------------------------*\

I added this section in the hope that it might spare quite some time to anyone
who will unfortunately have to deal with this code after me. It shows how the
1-D matrices in updateLocalTemperatureProfile are built.

Let us start with the heat diffusion equation in cylindrical coordinates with
temporally constant cp, rho and spatially constant k. For simplicity,
let us not consider a pin right now, just a regular uniform hollow cylinder
cylinder:

rho*cp*ddtT - div(k*grad(T)) = q    =>    rho*cp*ddt - (k/r)*ddr(r*ddr(T)) = q

Here is an example of a mesh with 5 nodes, represented by (o) and two "ghost"
nodes before and after the boundary nodes 0 and 4.

    radius=   rIn  ...                           ...  rOut

    index =   0         1         2         3         4 (= N-1)

    x_ _ _ _ _o_________o_________o_________o_________o_ _ _ _ _x



Let us discretise the eq. term by term with a finite difference scheme. For a 
radial node i, one has:

    1) rho*cp*ddtT -> rho_i*cp_i*( T_i - T_i_old )/dt =

        =       T_i *       ( rho_i*cp_i/dt )
            -   T_(i, old)* (rho_i*cp_i/dt)


    
    2) (k/r)*ddr(r*ddr(T)) = k*( (1/r)*(ddr(T)) + d2dr2(T) ) =>
        
        

        2.1) (k/r)*ddr(T) -> k*( T_(i+1) - T_(i-1) )/(2*dr*r_i) [CDS scheme] =

            =       T_(i+1) *   (  k/(2*dr*r_i) )
                +   T_(i-1) *   ( -k/(2*dr*r_i) )
        
        

        2.2) k*d2dr2(T) = ddr(ddr(T)) -> k*( ddr(T)_(i+1) - ddr(T)_i )/dr
            
            = k_i*( (T_(i+1) - T_(i))/dr - (T_(i) - T_(i-1))/dr )/dr =

            =       T_(i+1) *   (   k/(dr^2) )
                +   T_i *       ( -2k/(dr^2) )
                +   T_(i-1) *   (   k/(dr^2) )



        So (k/r)*ddr(r*ddr(T)) ->

            ->      T_(i+1) *   (   k*( 1/(dr^2)   + 1.0/(2*dr*r_i) ) )
                +   T_i *       ( -2k*( 1/(dr^2) )                    )
                +   T_(i-1) *   (   k*( 1/(dr^2)   - 1.0/(2*dr*r_i) ) )



    Recall that (k/r)*ddr(r*ddr(T)) appears with the - sign in the heat
    equation, so all the coefficients from 2) need to be changed in sign. 
    Combining the coefficients from 1) and 2) and 2.2) we get equation [I]:

        rho*cp*ddtT - div(k*grad(T)) = q ->

        ->      T_(i+1) *   ( -k*( 1/(dr^2)     + 1.0/(2*dr*r_i) ) )
            +   T_i *       ( 2k*( 1/(dr^2) )   + rho_i*cp_i/dt    )
            +   T_(i-1) *   ( -k*( 1/(dr^2)     - 1.0/(2*dr*r_i) ) )
            =
                q_i
            +   T_(i, old)*(rho_i*cp_i/dt)

    If we set:

        A = k/(dr^2)
        B = 1.0/(2*dr*r_i)
        X = rho_i*cp_i/dr

    Then we can re-write the whole thing as:                                 

        rho*cp*ddtT - div(k*grad(T)) = q ->

        ->      T_(i+1) *   ( -A + B )
            +   T_i *       ( 2A + X )
            +   T_(i-1) *   ( -A - B )
            =
                q_i
            +   T_(i, old)*X                                                [I]



Equation I is valid everywhere, even at the boundaries. However, at the
boundaries, T_N and T_(-1) do not exists. This is where ghost nodes come into
play. Let us have the following boundary conditions:

    1) Zero gradient at i = 0 : (i.e. Neumann)

        ddr(T)_(r_0) = 0 -> (T_1 - T(_(-1))/2*dr = 0 - > T_(-1) = T_1

    T(-1) does not exists and is our left ghost node. It does not matter,
    we obtained a relationship between T_1 and T_(-1)! Let us substitute this
    relationship in equation [I] and we get the equation at the left boundary:

        rho*cp*ddtT - div(k*grad(T)) = q -> (zeroGradient at r = rIn) ->

        ->      T_1 *       ( -2A )
            +   T_0 *       ( 2A + X )
            =
                q_i  
            +   T_(0, old)*X                                               [II]



    2) Convective boundary condition at i = N-1 : (i.e. Robin)

        k*ddr(T)_N = hExt*(TExt - T_(N-1))

    with hExt being the heat transfer coefficients of the medium after
    rOut and TExt being the temperature of the medium. As before, we can
    discretize it as:

        k*( T_N - T_(N-2) )/(2*dr) = hExt*(TExt - T_(N-1)) ->

        -> T_N = T_(N-2) + (2*dr/k)*hExt*(TExt -T_(N-1))

    If we substitute the expression for the ghost node T_N in equation I, we 
    obtain the BC. For convenience, let us rename C = (2*dr/k):

        rho*cp*ddtT - div(k*grad(T)) = q -> (convective at r = rOut) ->

                T_(N-1) * ( 2A + hExt*C*(A+B) + X )
            +   T_(N-2) * ( -2A )
            =
                q_(N-1)
            +   T_(N-1, old)*X
            +   TExt*hExt*C*(A+B)                                         [III]



So, there you have it, the boundary coefficients (II, III) and the bulk
coefficients (I).
Ok, now what if you want to do a cylinder? Well, you don't need to change
anything as there are no terms that contain 1/r_i for the zeroGradient
expression at r = r_1. Ok again, what if you want to do a proper nuclear fuel 
pin? Well, you should be able to understand what I did in the code. I 
used a convective boundary condition at the fuel outer surface where TExt is 
the  actual T_(i+1), while the fuel T_(i+1) is treated as a ghost node. Same 
for the inner side of the cladding, yet the coefficients are swapped for 
obvious reasons. For the outer cladding surface, again, it is a convective 
boundary conditions, yet it is complicated by the fact that I might have a mix 
of vapour and liquid, at different temperatures, both contacting the pin. For 
simplicity, I assume a single cladding outer temperature. Then, due to the
additivity of heat transfer phenomena, I simply consider:

    k*( T_N - T_(N-2) )/2*dr = 
        =   hVap*fracVap*(TVap - T_(N-1)) + hLiq*fracLiq*(TLiq - T_(N-1))

with fracVap and fracLiq being the fractions of structure interfacial area
that is in contact with either the vapour or the liquid. More in general, the
treatment of the convective boundary conditions can be generalized for any
number of fluids contacting the pin.
Let us have:

    k*( T_N - T_(N-2) )/2*dr = 
        =   SUM_j [hf_j(T_(ext, j) - T_(N-1))]

Where SUM_j denotes a summation over the j indices. hf_j is the heat
transfer coefficient betwee the j-th fluid and the structure, multiplied by
the fraction of interfacial area of the structure that is in contact with the
j-th fluid. Then, this can be re-written as:

    k*( T_N - T_(N-2) )/2*dr = 
        =   SUM_j[hf_j*T_(ext, j)] - SUM_j[hf_j]*T_(N-1)
        =   HTSum - HSum*T_(N-1)

which is the notation that is used in the code implementation. This grants more
generality to this class, which does not need to now any details on how many 
fluid are there, nor what are the interfacial area fractions and so on, as
these need to be passed by the user at a higher level (i.e., when calling the
correct(HTSum, HSum) function in the main program).

\*---------------------------------------------------------------------------*/