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

#include "twoPhaseDragMultiplierModel.H"
#include "fluidGeometry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseDragMultiplierModel, 0);
    defineRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel, 
        twoPhaseDragMultiplierModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragMultiplierModel::twoPhaseDragMultiplierModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("twoPhaseDragMultiplierModel", objReg.name()),
            mesh.time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    mFluid_
    (
        mesh_.lookupObject<fluid>
        (
            word("alpha."+dict.get<word>("multiplierFluid"))
        )
    ),
    mFracPtr_(nullptr),
    oFracPtr_(nullptr),
    phi2_
    (
        IOobject
        (
            IOobject::groupName("phi2", objReg.name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0.0)
    ),
    maxPhi2_
    (
        dict.lookupOrDefault<scalar>("maxValue", 1000)
    )
{
    HashTable<const fluid*> fluids(mesh_.lookupClass<fluid>());
    wordList fluidNames(fluids.toc());

    //- PORCO DIOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    oFluidPtr_ = (fluidNames[0] == "alpha."+mFluid_.name()) ?
        fluids[fluidNames[1]] : fluids[fluidNames[0]];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseDragMultiplierModel::setPhi2(const volScalarField& phi2)
{
    //- Under-relax if necessary
    scalar f(1.0);
    if (mesh_.relaxField("twoPhaseDragMultiplier"))
        f = mesh_.fieldRelaxationFactor("twoPhaseDragMultiplier");
    if (f!= 1.0) 
        phi2_ = f*phi2+(1.0-f)*phi2_;
    else
        phi2_ = phi2;

    //- Limit extrema
    phi2_ = max(min(phi2_, maxPhi2_), 1.0);  

    //- Limit spatially only to two-phase regions
    forAll(phi2_, i)
    {
        scalar phaseProduct(mFluid_[i]*(*oFluidPtr_)[i]);

        if (phaseProduct < 1e-6)
        {
            phi2_[i] = 1.0;
        }
    }  
}

/*
void Foam::twoPhaseDragMultiplierModel::limitPhi2
(
    const volScalarField& x,
    const scalar& x0,
    const scalar& x1
)
{
    this->limitPhi2();

    volScalarField marker(pos(x-x0));

    phi2_ = 
        marker*
        (
            phi2_
        +   (1.0-phi2_)*min((x-x0)/(x1-x0), 1.0)
        ) 
    +   (1.0-marker)*phi2_;
    
}
*/

/*

Here is the theoretical explanation. Suppose I have my two fluids, each having
a momentum conservation equation in the form:

ddt(alpha_i*rho_i*u_i) + ... = -F_wi + ...

With F_wi being the volumetric force due to structure drag acting on phase
i. In FFSEulerFoam:

F_wi = Kd_wi*u_i

with the general form of Kd_wi being:

Kd_wi = gamma*f_w(Re_i)*0.5*rho_i*mag(u_i)/Dh

with gamma being the structure void fraction (i.e. 1.0-alpha_s) and

Re_i = u_i*Dh/nu_i

being the Reynolds number.

The idea behind the twoPhaseFlow multiplier is saying that the TOTAL volumetric
force due to structure drag (i.e. the one arising when considering the two
phase mixture, i.e. the sum of momentum equations of both phases) is

F_wtot = phi2_I*F_wI1p

where phi2_ is the two phase multiplier and F_wI1p is the volumetric drag
force that would exist IF phase I would flow through the entirety of the
available structure cross-section (i.e. assuming that there is no other phase).
Please note that due to the porous context, even if there was no other phase,
phase I would still flow at most in a gamma fraction of the total available 
mesh volume. For simplicity, I will define alphaNorm_i = alpha_i/gamma. 
Ultimatley, the I1p subscript means "if "phase I flows throught the entirety of
the available void fraction as if it was a single-phase flow".
So, what does F_wI1p look like?

F_wI1p = Kd_wI1p*u_I1p

with u_I1p = alphaNorm_i*u_i (should be easy to see why, same reasoning as
the one used in the definition of the superficial velocity in porous media).

Kd_wI1p = gamma*f_wI(Re_I1p)*0.5*rho_I*mag(u_I1p)/Dh

With:

Re_I1p = u_i1p*Dh/nu_i = alphaNorm_i*u_i*Dh/nu_i

Please note that in purely single phase scenarios Re_I1p and Re_i are
exactly the same (alphaNorm_i = 1.0)! For this reason, all of the drag 
models are always computed with Re_I1p (so I get the right quantity in
both scenarios. Re is computed by the FSPair class, have a look there).

Ok, so let us go back to the discussion. I know how to compute the total
volumetric drag F_wtot starting from the F_wI1p drag of ONE of the two phases.
(which is the fluid specified by mFluid_ in this class).
How do I compute the individual F_wi components (whose sum is F_wtot) to
apply to the individual phase momentum equations? I use structure contact
fraction areas frac_i, which tell me which fluid is contacting the structure.
So:

F_wi = frac_i*F_wtot

with the frac_i being calculated by the structureInterfacialAreaPartitonModel
that applies in the regime of interest.

The final expression of F_wi thus is:

F_wi    =   frac_i*phi2_I*(gamma*f_w(Re_I1p)*0.5*rho_I*mag(u_I1p)/Dh)*u_I1p  ->
        
        =   frac_i*phi2_I*(alphaNorm_i)^2*
            (gamma*f_w(Re_I1p)*0.5*rho_I*mag(u_I)/Dh)*u_I 
        
        =   frac_i*phi2_I*(alphaNorm_i)^2*Kd_wI*u_I

With i != I! This is the crucial part! Let us make an example. Let us consider
two phases, 1 and 2. Let us assume that we want to compute F_w1, F_w2 from the
total volumetric F_wtot by using the phase 1 single phase pressure drop. We 
have:

F_w1 = frac_1*phi2_1*(alphaNorm_1)^2*Kd_w1*u_1
F_w2 = frac_2*phi2_1*(alphaNorm_1)^2*Kd_w1*u_1

We define modKd_wi:

modKd_wi = frac_i*phi2_I*(alphaNorm_I)^2*Kd_wI

So that:

F_w1 = modKd_1*u_1
F_w2 = modKd_2*u_1

The correctKdTable essential turns Kd_wi into modKd_wi. However, the drag is
finally added in the momentum equations (see UEqns.H) always and exclusively
as -modKd_i*u_i. This is fine for phase 1 (which is the one on which the
multiplier is based). This NOT fine for phase 2, as u_2 should be used, not
u_1. Thus, modKd_2 is additionally modified by mag(u_1)/mag(u_2). Now, is
this correct? Well, no, because I am losing directional information by 
rescaling by magnitudes, yet I am not bothered by this as the two phase
multiplier method was not even meant to be applied outside 1-D cases, yet
people still do it (see SABENA). To connect this example with the code 
below, KdI is Kd1 and KdO is Kd2.

*/

void Foam::twoPhaseDragMultiplierModel::correctKdTable
(
    const fluidGeometry& fG,
    volTensorFieldPtrTable& Kds
)
{
    //- Set ptrs to contact fractions if not set
    if (mFracPtr_ == nullptr)
    {
        mFracPtr_ = &(fG.fKd(mFluid_.name()));
    }
    if (oFracPtr_ == nullptr)
    {
        oFracPtr_ = &(fG.fKd(oFluidPtr_->name()));
    }
    const volScalarField& mFrac(*mFracPtr_);
    const volScalarField& oFrac(*oFracPtr_);

    volTensorField& KdI(*Kds[mFluid_.name()+".structure"]);
    volTensorField& KdO(*Kds[oFluidPtr_->name()+".structure"]);
    volScalarField corr(phi2_*sqr(mFluid_.normalized()));

    /*
    Info << mFluid_.name() << " " << oFluidPtr_->name() << endl;    
    Info << "Before" << endl;
    labelList cells(0);
    label i0(500);
    for (int i = i0; i<i0+20; i++)
    {
        cells.append(i);
    }
    forAll(cells, i)
    {
        label celli(cells[i]);
        Info << mFrac[i] << " " << KdI[celli] << " " << oFrac[i] << " " << KdO[celli] << " " << phi2_[celli] << endl;
    }
    */

    KdI *= corr;
    KdO =   
        KdI*mFluid_.magU()/
        (
            max
            (
                oFluidPtr_->magU(),
                dimensionedScalar("", dimVelocity, 1e-3)
            )
        );
    /*
    Info << "After1" << endl;
    forAll(cells, i)
    {
        label celli(cells[i]);
        Info << KdI[celli] << " " << KdO[celli] << endl;
    }
    */
    KdI *= mFrac;
    KdO *= oFrac;

    /*
    Info << "After2" << endl;
    forAll(cells, i)
    {
        label celli(cells[i]);
        Info << KdI[celli] << " " << KdO[celli] << endl;
    }
    */
}

// ************************************************************************* //
