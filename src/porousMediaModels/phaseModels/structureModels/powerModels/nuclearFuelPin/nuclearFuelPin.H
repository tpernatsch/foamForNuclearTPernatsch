/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::powerModels::nuclearFuelPin

Description
    Model for representing a nuclear fuel pin with constant material
    properties. Models fuel, gap, cladding and solves for the heat equation
    with a finite differences scheme. Only radial heat diffusion is considered.
    Derivation of the equations is provided at the end of the present file.

## Derivation

\warning
THIS IS THE DERIVATION OF THE PREVIOUS FINITE-DIFFERENCE EQUATIONS
WHICH ARE NO LONGER USED. THE CURRENT CODE IS BASED ON FINITE VOLUME SCHEME
WHICH IS SIGNIFICANTLY MORE CONSERVATIVE THAN THE FINITE-DIFFERENCE ONE
FOR THE SAME MESH SIZE. THE FINITE VOLUME DERIVATION IS NOT REPORTED AS...

I added this section in the hope that it might spare quite some time to anyone
who will unfortunately have to deal with this code after me. It shows how the
1-D matrices in updateLocalTemperatureProfile are built.

Let us start with the heat diffusion equation in cylindrical coordinates with
temporally constant cp, rho and spatially constant k. For simplicity,
let us not consider a pin right now, just a regular uniform hollow cylinder:

$$
\rho c_p \frac{dT}{dt} - div(k grad(T)) = q    \rightarrow    \rho c_p \frac{dT}{dt} - \frac{k}{r} \frac{d}{dr} r \frac{dT}{dr}  = q
$$

Here is an example of a mesh with 5 nodes, represented by (o) and two "ghost"
nodes before and after the boundary nodes 0 and 4.

```
    radius=   rIn  ...                           ...  rOut

    index =   0         1         2         3         4 (= N-1)

    x_ _ _ _ _o_________o_________o_________o_________o_ _ _ _ _x
```


Let us discretise the eq. term by term with a finite difference scheme. For a
radial node $i$, one has:

1)

$$
\rho c_p \frac{dT}{dt} \rightarrow
    \rho_i c_{p_i} \frac{T_i - T_{i,\text{old}}}{dt} =
    T_i \left( \frac{\rho_i c_{p,i}}{dt} \right)
    - T_{i,\text{old}} \left( \frac{\rho_i c_{p,i}}{dt} \right)
$$

2)

$$
\frac{k}{r} \frac{d}{dr} r \frac{dT}{dr} = k \left( \frac{1}{r}\frac{dT}{dr} + \frac{d^2T}{dr^2} \right)
$$

2.1)

$$
\frac{k}{r}\frac{dT}{dr} \rightarrow k \frac{T_{i+1} - T_{i-1}}{2 r_i dr} \text{[CDS scheme]} =
      T_{i+1} \frac{k}{2 r_i dr}
    - T_{i-1} \frac{k}{2 r_i dr}
$$

2.2)

$$
\begin{align}
k \frac{d^2T}{dr^2} = \frac{d}{dr} \frac{dT}{dr} \rightarrow
k \frac{ \frac{dT_{i+1} }{dr} - \frac{dT_i}{dr} }{dr}
    & = k_i \frac{ \frac{T_{i+1} - T_i}{dr} - \frac{T_i - T_{i-1}}{dr} }{dr} \\
    & = T_{i+1} \frac{k}{dr^2}
        - T_i \frac{2k}{dr^2}
        + T_{i-1} \frac{k}{dr^2}
\end{align}
$$

So

$$
\frac{k}{r} \frac{d}{dr} r \frac{dT}{dr} \rightarrow
    T_{i+1} k \left( \frac{1}{dr^2} + \frac{1}{2 r_i dr} \right)
-   T_i \frac{2k}{dr^2}
+   T_{i-1} k \left( \frac{1}{dr^2} - \frac{1}{2 r_i dr} \right)
$$


Recall that $\frac{k}{r} \frac{d}{dr} r \frac{dT}{dr}$ appears with the - sign in the heat
equation, so all the coefficients from 2) need to be changed in sign.
Combining the coefficients from 1) and 2) and 2.2) we get equation [I]:

$$
\begin{align}
& \rho c_p \frac{dT}{dt} - div(k grad(T)) = q \\
& \rightarrow
-   T_{i+1} k \left( \frac{1}{dr^2}  - \frac{1}{2 r_i dr}        \right)
+   T_i       \left( \frac{2k}{dr^2} + \frac{\rho_i c_{p,i}}{dt} \right)
-   T_{i-1} k \left( \frac{1}{dr^2}  + \frac{1}{2 r_i dr}        \right)
=
    q_i
+   T_{i, old} \frac{\rho_i c_{p,i}}{dt}
\end{align}
$$

If we set:

$$
\begin{align}
A & = \frac{k}{dr^2} \\
B & = \frac{1}{2 r_i dr} \\
X & = \frac{\rho_i c_{p,i}}{dt}
\end{align}
$$

Then we can re-write the whole thing as (referred later as eq. I):

$$
\begin{align}
& \rho c_p \frac{dT}{dt} - div(k grad(T)) = q \\
& \rightarrow
    T_{i+1} ( -A + B )
+   T_i     ( 2A + X )
+   T_{i-1} ( -A - B )
=
    q_i
+   T_{i, old} X
\end{align}
$$


Equation I is valid everywhere, even at the boundaries. However, at the
boundaries, $T_N$ and $T_{-1}$ do not exists. This is where ghost nodes come into
play. Let us have the following boundary conditions:

    1) Zero gradient at $i = 0$ : (i.e. Neumann)

$$
\begin{align}
\frac{dT}{dr} |_{r_0} = 0  & \rightarrow  \frac{T_1 - T_{-1}}{2} dr = 0 \\
                           & \rightarrow  T_{-1} = T_1
\end{align}
$$

    $T_{-1}$ does not exists and is our left ghost node. It does not matter,
    we obtained a relationship between $T_1$ and $T_{-1}$! Let us substitute this
    relationship in equation [I] and we get the equation at the left boundary
    (referred later as eq. [II]):

$$
\begin{align}
\rho c_p \frac{dT}{dt} - div(k grad(T)) = q
    & \rightarrow (\text{zeroGradient at} r = r_{In}) \\
    & \rightarrow T_1 \times ( -2A ) + T_0 \times ( 2A + X ) = q_i + T_{0, old} X
\end{align}
$$

    2) Convective boundary condition at $i = N-1$ : (i.e. Robin)

$$
k \frac{dT}{dr} |_N = h_{Ext} (T_{Ext} - T_{N-1})
$$

    with hExt being the heat transfer coefficients of the medium after
    rOut and TExt being the temperature of the medium. As before, we can
    discretize it as:

$$
k \frac{ T_N - T_{N-2} }{2 dr} = h_{Ext} (T_{Ext} - T_{N-1})
\rightarrow
T_N = T_{N-2} + \frac{2 dr}{k} h_{Ext} (T_{Ext} - T_{N-1})
$$

    If we substitute the expression for the ghost node T_N in equation I, we
    obtain the BC. For convenience, let us rename $C = \frac{2 dr}{k}$
    (referred later as eq. [III]):

$$
\begin{align}
\rho c_p \frac{dT}{dt} - div(k grad(T)) = q
    & \rightarrow (\text{convective at} r = r_{Out}) \\
    & \rightarrow
                T_{N-1} ( 2A + h_{Ext} C (A+B) + X )
            +   T_{N-2} \times ( -2A )
            =
                q_{N-1}
            +   T_{N-1, old} X
            +   T_{Ext} h_{Ext} C (A+B)
\end{align}
$$

So, there you have it, the boundary coefficients (II, III) and the bulk
coefficients (I).
Ok, now what if you want to do a cylinder? Well, you don't need to change
anything as there are no terms that contain 1/r_i for the zeroGradient
expression at $r = r_1$. Ok again, what if you want to do a proper nuclear fuel
pin? Well, you should be able to understand what I did in the code. I
used a convective boundary condition at the fuel outer surface where TExt is
the  actual $T_{i+1}$, while the fuel $T_{i+1}$ is treated as a ghost node. Same
for the inner side of the cladding, yet the coefficients are swapped for
obvious reasons. For the outer cladding surface, again, it is a convective
boundary conditions, yet it is complicated by the fact that I might have a mix
of vapour and liquid, at different temperatures, both contacting the pin. For
simplicity, I assume a single cladding outer temperature. Then, due to the
additivity of heat transfer phenomena, I simply consider:

$$
k \frac{T_N - T_{N-2}}{2} dr =
    h_{Vap} f_{Vap} (T_{Vap} - T_{N-1}) + h_{Liq} f_{Liq} (T_{Liq} - T_{N-1})
$$

with $f_{Vap}$ and $f_{Liq}$ being the fractions of structure interfacial area
that is in contact with either the vapour or the liquid. More in general, the
treatment of the convective boundary conditions can be generalized for any
number of fluids contacting the pin.
Let us have:

$$
k \frac{T_N - T_{N-2}}{2} dr = \sum_j h_{f,j} (T_{ext, j} - T_{N-1})
$$

$h_{f,j}$ is the heat
transfer coefficient betwee the j-th fluid and the structure, multiplied by
the fraction of interfacial area of the structure that is in contact with the
j-th fluid. Then, this can be re-written as:

$$
\begin{align}
k \frac{T_N - T_{N-2}}{2} dr & = \sum_j h_{f,j} T_{ext, j} - T_{N-1} \sum_j h_{f,j} \\
                             & = HTSum - T_{N-1} \times HSum
\end{align}
$$

which is the notation that is used in the code implementation. This grants more
generality to this class, which does not need to now any details on how many
fluid are there, nor what are the interfacial area fractions and so on, as
these need to be passed by the user at a higher level (i.e., when calling the
correct(HTSum, HSum) function in the main program).

\par Options

    Parameters in the powerModel subdictionary for `nuclearFuelPin`:

    - <b>`powerDensity`</b> - Initial uniform power density if not computed by
    the neutronics sub-solver.
    - <b>`fuelInnerRadius`</b> - Fuel pellet inner radius.
    - <b>`fuelOuterRadius`</b> - Fuel pellet outer radius.
    - <b>`cladInnerRadius`</b> - Cladding inner radius.
    - <b>`cladOuterRadius`</b> - Cladding outer radius.
    - <b>`fuelMeshSize`</b> - Fuel pellet radial mesh size.
    - <b>`cladMeshSize`</b> - Cladding radial mesh size.
    - <b>`fuelRho`</b> - Fuel constant density in kg/m3.
    - <b>`fuelCp`</b> - Fuel constant specific heat capacity in J/kg/K.
    - <b>`fuelK`</b> - Fuel constant thermal conductivity in W/m/K.
    - <b>`cladRho`</b> - Cladding constant density in kg/m3.
    - <b>`cladCp`</b> - Cladding constant specific heat capacity in J/kg/K.
    - <b>`cladK`</b> - Cladding constant thermal conductivity in W/m/K.
    - <b>`fuelT`</b> - Fuel initial temperature.
    - <b>`cladT`</b> - Cladding initial temperature.
    - <b>`gapH`</b> - Gap gas constant heat transfer coefficient.
    - <b>`gapHPowerDensity`</b> - Gap gas power-dependent heat transfer
    coefficient. Can be used in replacement of `gapH`.

Usage
    From the one phase and PK tutorial.
    \verbatim
        powerModel
        {
            type                nuclearFuelPin;

            // Initial power density if not computed by the neutronics sub-solver
            powerDensity        0;

            // Fuel pin dimensions in m
            fuelInnerRadius     0.0012;
            fuelOuterRadius     0.004715;
            cladInnerRadius     0.004865;
            cladOuterRadius     0.005365;

            // Fuel pin radial mesh
            fuelMeshSize        16;
            cladMeshSize        4;

            // Fuel constant thermophysical properties
            fuelRho             10480; // kg/m3
            fuelCp              250; // J/kg/K
            fuelK               3; // W/m/K

            // Cladding constant thermophysical properties
            cladRho             7500;
            cladCp              500;
            cladK               20;

            // Gap gas heat transfer coefficient
            gapH                3000;

            // Initial fuel and cladding temperatures
            fuelT               600;
            cladT               600;
        }
    \endverbatim

    From the FFTF tutorial.
    \verbatim
        powerModel
        {
            //- This volume fraction accounts for the pins only
            type                nuclearFuelPin;
            volumeFraction      0.5539729557;
            fuelMeshSize        10;
            cladMeshSize        5;
            fuelT               574.56;
            cladT               574.56;
            cladK               20;
            fuelK               3;
            cladRho             7870;
            fuelRho             10480;
            cladCp              500;
            fuelCp              250;
            cladOuterRadius     0.002921;
            cladInnerRadius     0.002540;
            fuelOuterRadius     0.002470;
            fuelInnerRadius     0.0;

            // Optional power-dependent heat transfer coefficient
            gapHPowerDensity    table
            (
                (   0               528.4782467867  )
                (   795443280.75916 3170.8694807204 )
            );
        }
    \endverbatim

SourceFiles
    nuclearFuelPin.C

\mainauthor
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com),
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)

\*---------------------------------------------------------------------------*/

#ifndef nuclearFuelPin_H
#define nuclearFuelPin_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class nuclearFuelPin Declaration
\*---------------------------------------------------------------------------*/

class nuclearFuelPin
:
    public powerModel
{
protected:

    //- Field (over the global mesh) of scalarFields (over a 1-D mesh of size
    //  subMeshSize_) representing the radial fuel temperature profile
    //  across fuel and cladding
    IOFieldField<Field, scalar> Trad_;

    //- Field representing pin centerlinetemperature
    volScalarField Tfi_;

    //- Field representing pin outer fuel temperature
    volScalarField Tfo_;

    //- Field representing pin inner cladding temperature
    volScalarField Tci_;

    //- Field representing pin outer cladding temperature
    volScalarField Tco_;

    //- Average fuel temperature
    volScalarField Tfav_;

    //- Average cladding temperature
    volScalarField Tcav_;

    //- Scalar that are input in the IOdictionary for passing max fuel
    //  temperature
    scalar Tfmax_;

    //- Scalar that are input in the IOdictionary for passing min fuel
    //  temperature
    scalar Tfmin_;

    //- Scalar that are input in the IOdictionary for passing max cladding
    //  temperature
    scalar Tcmax_;

    //- Scalar that are input in the IOdictionary for passing min cladding
    //  temperature
    scalar Tcmin_;

    //- Fraction of power from neutronics that goes to this powerModel
    //  This take sinto account for instance cases where some of the
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- Number of cells in the 1-D radial fuel mesh
    labelList fuelMeshSize_;

    //- Number of cells in the 1-D radial clad mesh
    labelList cladMeshSize_;

    //- Sum of fuelMeshSize and cladMeshSize
    labelList meshSize_;

    //- Radial position of the mesh nodes
    scalarListList r_;

    //- Fuel inner radius
    scalarList rfi_;

    //- Fuel outer radius
    scalarList rfo_;

    //- Cladding inner radius
    scalarList rci_;

    //- Cladding outer radius
    scalarList rco_;

    //- Delta radius in fuel (i.e. mesh cell length)
    scalarList drf_;

    //- Delta radius in cladding (i.e. mesh cell length)
    scalarList drc_;

    //- Gap width
    scalarList drg_;

    //- Fuel volumetric heat capacity
    scalarList rhoCpf_;

    //- Cladding volumetric heat capacity
    scalarList rhoCpc_;

    //- Fuel thermal conductivity
    scalarList kf_;

    //- Cladding thermal conductivity
    scalarList kc_;

    //- Gap heat transfer coefficient
    scalarList gapH_;

    //- True if fuel has hole
    boolList hollowFuel_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;

    //-
    wordList regionNames_;

    //- Just to avoid typing constant::mathematical::pi every time...
    scalar pi_;

    //- Mesh ring areas
    scalarListList dA_;

    //- Optional table to describe the dependence of the gap condunctance on
    //  the volumetric power density
    HashTable<autoPtr<Function1<scalar>>, word, word::hash>
        gapHPowerDensityTable_;

    //-
    boolList useGapHPowerDensityTable_;

    // - List to fuel fractions (i.e., fuel only, without considering the cladding)
    scalarList fuelFraction_;

public:

    //- Runtime type information
    TypeName("nuclearFuelPin");

    // Constructors

        nuclearFuelPin
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~nuclearFuelPin();

    // Member Functions

        //- Set interfacial area consistently with the pin geometry. In
        //  particular, if the overall total volume-fraction of the pins
        //  in a given mesh cell is alpha, it can be shown that the interfacial
        //  area density iA of the pin surface is 2*alpha/outerRadius
        virtual void setInterfacialArea();

        void updateLocalAvgGlobalMinMaxT
        (
            const label& starti,
            const label& endi,
            const scalarList& r,
            const scalar& dr,
            const scalarField& Trad,
            scalar& Tav,
            scalar& Tmin,
            scalar& Tmax
        );

        //- Update temperature profiles Tc, Tf, in a single domain cell with
        //  current system conditions
        void updateLocalTemperatureProfile
        (
            const label& celli,
            const scalar& ht,
            const scalar& h
        );

        //- In two-phase solvers :
        //  HT = H1*frac1*Tfluid1 + H2*frac2*Tfluid2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*Tfluid
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        );

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual void correctT(volScalarField& T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoSubstructures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //