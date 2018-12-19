/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "constSubscaleFuel.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"
#include "Field.H"
#include "scalarField.H"
#include "SquareMatrix.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constSubscaleFuel, 0);

    addToRunTimeSelectionTable
    (
        subscaleFuel,
        constSubscaleFuel,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constSubscaleFuel::constSubscaleFuel
(
    const fvMesh& mesh
)
:
    subscaleFuel(mesh),
    dFuel_(fuelZoneNumber_),
    dClad_(fuelZoneNumber_),
    planar_(fuelZoneNumber_),
    fuelVolPower_(fuelZoneNumber_),
    claddingK_(fuelZoneNumber_),
    fuelK_(fuelZoneNumber_),
    claddingRho_(fuelZoneNumber_),
    fuelRho_(fuelZoneNumber_),
    claddingCp_(fuelZoneNumber_),
    fuelCp_(fuelZoneNumber_),
    gapH_(fuelZoneNumber_),
    rcOut_(fuelZoneNumber_),
    rcIn_(fuelZoneNumber_),
    rfOut_(fuelZoneNumber_),
    rfIn_(fuelZoneNumber_)
{

    PtrList<entry> entries(IOdictionary::lookup("zones"));
    forAll(entries,zoneI)
    {

        dictionary& dict = entries[zoneI].dict();
        fuelVolPower_.set(zoneI,new scalar(dict.lookupOrDefault("fuelVolPower",0.0)));
        claddingK_.set(zoneI,new scalar(dict.lookupOrDefault("claddingK",0.0)));
        fuelK_.set(zoneI,new scalar(dict.lookupOrDefault("fuelK",0.0)));
        claddingRho_.set(zoneI,new scalar(dict.lookupOrDefault("claddingRho",0.0)));
        fuelRho_.set(zoneI,new scalar(dict.lookupOrDefault("fuelRho",0.0)));
        claddingCp_.set(zoneI,new scalar(dict.lookupOrDefault("claddingCp",0.0)));
        fuelCp_.set(zoneI,new scalar(dict.lookupOrDefault("fuelCp",0.0)));
        gapH_.set(zoneI,new scalar(dict.lookupOrDefault("gapH",0.0)));
        rcOut_.set(zoneI,new scalar(dict.lookupOrDefault("rcOut",0.0)));
        rcIn_.set(zoneI,new scalar(dict.lookupOrDefault("rcIn",0.0)));
        rfOut_.set(zoneI,new scalar(dict.lookupOrDefault("rfOut",0.0)));
        rfIn_.set(zoneI,new scalar(dict.lookupOrDefault("rfIn",0.0)));
        dFuel_.set(zoneI,new scalar((rfOut_[zoneI]-rfIn_[zoneI])/(fuelSubMeshSize_[zoneI]-1)));
        dClad_.set(zoneI,new scalar((rcOut_[zoneI]-rcIn_[zoneI])/(cladSubMeshSize_[zoneI]-1)));
        planar_.set(zoneI,new bool(dict.lookupOrDefault("planar",false)));

    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constSubscaleFuel::~constSubscaleFuel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::constSubscaleFuel::updateLocalHeatFluxImplicit(Foam::label zoneI, Foam::label cellIlocal, Foam::scalar Tfl, Foam::scalar hfl, Foam::scalar q)
{
    Field<scalar>& TcOld = Tc_[zoneI].oldTime()[cellIlocal];
    Field<scalar>& TfOld = Tf_[zoneI].oldTime()[cellIlocal];

    scalar dC = dClad_[zoneI];
    scalar kC = claddingK_[zoneI];
    scalar alphaClad = (kC/(claddingRho_[zoneI]*claddingCp_[zoneI]));
    scalar rC = alphaClad * mesh_.time().deltaT().value() /  pow(dC,2);
    scalar BetaMC = 1 + (1 + 1/(2*(rcOut_[zoneI]/dC))) * dC * hfl / kC; //(cladSubMeshSize_[zoneI]-1+rcIn_[zoneI]/dC)
    scalar gammaMC = (1 + 1/(2*(rcOut_[zoneI]/dC))) * dC * hfl * Tfl / kC; //(cladSubMeshSize_[zoneI]-1+rcIn_[zoneI]/dC)
    scalar Beta0C = 1 + (1 - 1/(2*(rcIn_[zoneI]/dC))) * dC * gapH_[zoneI] * (rfOut_[zoneI]/rcIn_[zoneI]) / kC;
    scalar gamma0C = (1 - 1/(2*(rcIn_[zoneI]/dC))) * dC * gapH_[zoneI] * (rfOut_[zoneI]/rcIn_[zoneI]) * TfOld[fuelSubMeshSize_[zoneI]-1] / kC;

    scalar dF = dFuel_[zoneI];
    scalar kF = fuelK_[zoneI];
    scalar alphaFuel = (kF/(fuelRho_[zoneI]*fuelCp_[zoneI]));
    scalar rF = alphaFuel * mesh_.time().deltaT().value() /  pow(dF,2);
    scalar BetaMF = 1 + (1 + 1/(2*(rfOut_[zoneI]/dF))) * dF * gapH_[zoneI] / kF; //(fuelSubMeshSize_[zoneI]-1+rfIn_[zoneI]/dF)
    scalar gammaMF = (1 + 1/(2*(rfOut_[zoneI]/dF))) * dF * gapH_[zoneI] * TcOld[0] / kF; //(fuelSubMeshSize_[zoneI]-1+rfIn_[zoneI]/dF)
    scalar GF = pow(dF,2) * (q + fuelVolPower_[zoneI]) / kF ;
    scalar Beta0F = 1 ;

    SquareMatrix<scalar> fuelMatrix(fuelSubMeshSize_[zoneI],fuelSubMeshSize_[zoneI],Foam::zero());
    List<scalar> fuelSource(fuelSubMeshSize_[zoneI],0.0);
    SquareMatrix<scalar> cladMatrix(cladSubMeshSize_[zoneI],cladSubMeshSize_[zoneI],Foam::zero());
    List<scalar> cladSource(cladSubMeshSize_[zoneI],0.0);

    forAll(cladSource,i)
    {

        if(i == 0)
        {
                cladMatrix[i][i] = -2 * rC * Beta0C - 1;
                cladMatrix[i][i+1] = 2 * rC;

                cladSource[i] =  -TcOld[i] - 2 * rC * gamma0C ;
        }
        else if(i == (cladSubMeshSize_[zoneI] - 1))
        {
            {
                cladMatrix[i][i-1] = 2 * rC;
                cladMatrix[i][i] = -2 * rC * BetaMC - 1;

                cladSource[i] =  -TcOld[i] - 2 * rC * gammaMC ;
            }
        }
        else
        {
                cladMatrix[i][i-1] = rC * (1.0 - 1/(2*(i+rcIn_[zoneI]/dC)));
                cladMatrix[i][i] = -2 * rC - 1;
                cladMatrix[i][i+1] = rC * (1.0 + 1/(2*(i+rcIn_[zoneI]/dC)));

                cladSource[i] =  -TcOld[i] ;
        }
    }

    forAll(fuelSource,i)
    {

        if(i == 0)
        {
            if(rfIn_[zoneI] == 0.0)
            {
                fuelMatrix[i][i] = -4 * rF - 1;
                fuelMatrix[i][i+1] = 4 * rF;

                fuelSource[i] =  -TfOld[i] - rF * GF ;
            }
            else
            {
                fuelMatrix[i][i] = -2 * rF * Beta0F  - 1;
                fuelMatrix[i][i+1] = 2 * rF;

                fuelSource[i] =  -TfOld[i] - rF * GF ;
            }
        }
        else if(i == (fuelSubMeshSize_[zoneI] - 1))
        {
            fuelMatrix[i][i-1] = 2 * rF;
            fuelMatrix[i][i] = -2 * rF * BetaMF  - 1;

            fuelSource[i] =  -TfOld[i] - rF * GF - 2 * rF * gammaMF ;
        }
        else
        {
            fuelMatrix[i][i-1] = rF * (1.0 - 1/(2*(i+rfIn_[zoneI]/dF)));
            fuelMatrix[i][i] = -2 * rF - 1;
            fuelMatrix[i][i+1] = rF * (1.0 + 1/(2*(i+rfIn_[zoneI]/dF)));

            fuelSource[i] =  -TfOld[i] - rF * GF ;
        }
    }

    solve(Tc_[zoneI][cellIlocal],cladMatrix,cladSource);
    solve(Tf_[zoneI][cellIlocal],fuelMatrix,fuelSource);

    scalar qOut
    (
        hfl*
        (
            Tc_[zoneI][cellIlocal][cladSubMeshSize_[zoneI]-1] -
            Tfl
        )
    );

    
    return qOut;
}



Foam::scalar
Foam::constSubscaleFuel::updateLocalHeatFluxImplicitPlanar(Foam::label zoneI, Foam::label cellIlocal, Foam::scalar Tfl, Foam::scalar hfl, Foam::scalar q)
{
    Field<scalar>& TcOld = Tc_[zoneI].oldTime()[cellIlocal];
    Field<scalar>& TfOld = Tf_[zoneI].oldTime()[cellIlocal];

    scalar dC = dClad_[zoneI];
    scalar kC = claddingK_[zoneI];
    scalar alphaClad = (kC/(claddingRho_[zoneI]*claddingCp_[zoneI]));
    scalar rC = alphaClad * mesh_.time().deltaT().value() /  pow(dC,2);
    scalar BetaMC = 1 + dC * hfl / kC;
    scalar gammaMC = dC * hfl * Tfl / kC;
    scalar Beta0C = 1 + dC * gapH_[zoneI] / kC;
    scalar gamma0C = dC * gapH_[zoneI] * TfOld[fuelSubMeshSize_[zoneI]-1] / kC;

    scalar dF = dFuel_[zoneI];
    scalar kF = fuelK_[zoneI];
    scalar alphaFuel = (kF/(fuelRho_[zoneI]*fuelCp_[zoneI]));
    scalar rF = alphaFuel * mesh_.time().deltaT().value() /  pow(dF,2);
    scalar BetaMF = 1 + dF * gapH_[zoneI] / kF;
    scalar gammaMF = dF * gapH_[zoneI] * TcOld[0] / kF;
    scalar GF = pow(dF,2) * (q + fuelVolPower_[zoneI]) / kF ;
    scalar Beta0F = 1 ;

    SquareMatrix<scalar> fuelMatrix(fuelSubMeshSize_[zoneI],fuelSubMeshSize_[zoneI],Foam::zero());
    List<scalar> fuelSource(fuelSubMeshSize_[zoneI],0.0);
    SquareMatrix<scalar> cladMatrix(cladSubMeshSize_[zoneI],cladSubMeshSize_[zoneI],Foam::zero());
    List<scalar> cladSource(cladSubMeshSize_[zoneI],0.0);

    forAll(cladSource,i)
    {

        if(i == 0)
        {
            cladMatrix[i][i] = -2 * rC * Beta0C - 1;
            cladMatrix[i][i+1] = 2 * rC;

            cladSource[i] =  -TcOld[i] - 2 * rC * gamma0C ;

        }
        else if(i == (cladSubMeshSize_[zoneI] - 1))
        {
            cladMatrix[i][i-1] = 2 * rC;
            cladMatrix[i][i] = -2 * rC * BetaMC - 1;

            cladSource[i] =  -TcOld[i] - 2 * rC * gammaMC ;
        }
        else
        {
            cladMatrix[i][i-1] = rC ;
            cladMatrix[i][i] = -2 * rC - 1;
            cladMatrix[i][i+1] = rC ;

            cladSource[i] =  -TcOld[i] ;
        }
    }

    forAll(fuelSource,i)
    {

        if(i == 0)
        {
            fuelMatrix[i][i] = -4 * rF - 1;
            fuelMatrix[i][i+1] = 4 * rF;

            fuelSource[i] =  -TfOld[i] - rF * GF ;

        }
        else if(i == (fuelSubMeshSize_[zoneI] - 1))
        {

            fuelMatrix[i][i-1] = 2 * rF;
            fuelMatrix[i][i] = -2 * rF * BetaMF  - 1;

            fuelSource[i] =  -TfOld[i] - rF * GF - 2 * rF * gammaMF ;

        }
        else
        {
            fuelMatrix[i][i-1] = rF ;
            fuelMatrix[i][i] = -2 * rF - 1;
            fuelMatrix[i][i+1] = rF ;

            fuelSource[i] =  -TfOld[i] - rF * GF ;
        }
    }

    solve(Tc_[zoneI][cellIlocal],cladMatrix,cladSource);
    solve(Tf_[zoneI][cellIlocal],fuelMatrix,fuelSource);

    return hfl * (Tc_[zoneI][cellIlocal][cladSubMeshSize_[zoneI]-1] - Tfl) ;
}

// Solve for fuel heat sources in fluid given the fluid temperature (Tfl), the heat transfer coefficient between fluid and fuel, and heat source in fuel (qs)
// qs is optional. One can use a constant heat source specified in the dictionary
Foam::tmp<Foam::volScalarField>
Foam::constSubscaleFuel::heatSources(const volScalarField& Tfl, const volScalarField& hfl, const volScalarField& qs)
{

    // Field to return
    tmp<volScalarField> tqf(
        new volScalarField
        (
            IOobject
            (
                "tqf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(1,0,-3,0,0,0,0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& qf = tqf.ref();

    // Calculate fuel anc clad temperatures and local heat flux on surface of clad
    forAll(Tf_,zoneI)
    {
        forAll(Tf_[zoneI], cellIlocal)
        {
            label zoneId(fuelZoneID_[zoneI]);

            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id

            // Solve system of equations to get fuel temperatures
            if(planar_[zoneI])
            {
                qf[cellIglobal] = updateLocalHeatFluxImplicitPlanar(zoneI,cellIlocal,Tfl[cellIglobal], hfl[cellIglobal], qs[cellIglobal]); // Pass local temperatures and solve for heat flux on fuel surface
            }
            else
            {
                qf[cellIglobal] = updateLocalHeatFluxImplicit(zoneI,cellIlocal,Tfl[cellIglobal], hfl[cellIglobal], qs[cellIglobal]); // Pass local temperatures and solve for heat flux on fuel surface
            }
        }

    }

    return tqf;
}


Foam::tmp<Foam::volScalarField>
Foam::constSubscaleFuel::fuelIn()
{
    tmp<volScalarField> tTinFuel(
        new volScalarField
        (
            IOobject
            (
                "TinFuel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& TinFuel = tTinFuel.ref();
    forAll(Tf_,zoneI)
    {
        forAll(Tf_[zoneI], cellIlocal)
        {
            label zoneId(fuelZoneID_[zoneI]);
            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
            TinFuel[cellIglobal] = Tf_[zoneI][cellIlocal][0];
        }
    }
    TinFuel.correctBoundaryConditions();
    return tTinFuel;
}

void
Foam::constSubscaleFuel::fuelAverage(volScalarField& TavFuel)
{
    forAll(Tf_,zoneI)
    {
        if(planar_[zoneI])
        {
            scalar dF = dFuel_[zoneI];
            forAll(Tf_[zoneI], cellIlocal)
            {
                label zoneId(fuelZoneID_[zoneI]);
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
                TavFuel[cellIglobal] = 0.0;
    
                forAll(Tf_[zoneI][cellIlocal],i)
                {
                    if(i == 0 || i == (fuelSubMeshSize_[zoneI] - 1))
                    {
                        TavFuel[cellIglobal] += Tf_[zoneI][cellIlocal][i]  * dF /2.0 ;
                    }
                    else
                    {
                        TavFuel[cellIglobal] += Tf_[zoneI][cellIlocal][i] * dF;
                    }
                }
                TavFuel[cellIglobal] *= 1 / (rfOut_[zoneI] - rfIn_[zoneI]);
            }
        }
        else
        {
            scalar dF = dFuel_[zoneI];
            forAll(Tf_[zoneI], cellIlocal)
            {
                label zoneId(fuelZoneID_[zoneI]);
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
                TavFuel[cellIglobal] = 0.0;
    
                forAll(Tf_[zoneI][cellIlocal],i)
                {
                    if(i == 0 || i == (fuelSubMeshSize_[zoneI] - 1))
                    {
                        TavFuel[cellIglobal] += Tf_[zoneI][cellIlocal][i] * (dF * i + rfIn_[zoneI]) * dF /2.0 ;
                    }
                    else
                    {
                        TavFuel[cellIglobal] += Tf_[zoneI][cellIlocal][i] * (dF * i + rfIn_[zoneI]) * dF;
                    }
                }
                TavFuel[cellIglobal] *= 2 / (pow(rfOut_[zoneI],2) - pow(rfIn_[zoneI],2));
            }
        }
    }
    TavFuel.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField>
Foam::constSubscaleFuel::fuelOut()
{
    tmp<volScalarField> tToutFuel(
        new volScalarField
        (
            IOobject
            (
                "ToutFuel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& ToutFuel = tToutFuel.ref();
    forAll(Tf_,zoneI)
    {
        forAll(Tf_[zoneI], cellIlocal)
        {
            label zoneId(fuelZoneID_[zoneI]);
            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
            ToutFuel[cellIglobal] = Tf_[zoneI][cellIlocal][fuelSubMeshSize_[zoneI]-1];
        }
    }
    ToutFuel.correctBoundaryConditions();
    return tToutFuel;
}


Foam::tmp<Foam::volScalarField>
Foam::constSubscaleFuel::cladIn()
{
    tmp<volScalarField> tTinClad(
        new volScalarField
        (
            IOobject
            (
                "TinClad",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& TinClad = tTinClad.ref();
    forAll(Tc_,zoneI)
    {
        forAll(Tc_[zoneI], cellIlocal)
        {
            label zoneId(fuelZoneID_[zoneI]);
            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
            TinClad[cellIglobal] = Tc_[zoneI][cellIlocal][0];
        }
    }
    TinClad.correctBoundaryConditions();
    return tTinClad;
}

void
Foam::constSubscaleFuel::cladAverage(volScalarField& TavClad)
{
    forAll(Tc_,zoneI)
    {
        if(planar_[zoneI])
        {
            scalar dC = dClad_[zoneI];
            forAll(Tc_[zoneI], cellIlocal)
            {
                label zoneId(fuelZoneID_[zoneI]);
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
                TavClad[cellIglobal] = 0.0;

                forAll(Tc_[zoneI][cellIlocal],i)
                {
                    if(i == 0 || i == (cladSubMeshSize_[zoneI] - 1))
                    {
                        TavClad[cellIglobal] += Tc_[zoneI][cellIlocal][i] * dC /2.0;
                    }
                    else
                    {
                        TavClad[cellIglobal] += Tc_[zoneI][cellIlocal][i] * dC;
                    }
                }
                TavClad[cellIglobal] *= 1 / (rcOut_[zoneI] - rcIn_[zoneI]);
            }
        }
        else
        {
            scalar dC = dClad_[zoneI];
            forAll(Tc_[zoneI], cellIlocal)
            {
                label zoneId(fuelZoneID_[zoneI]);
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
                TavClad[cellIglobal] = 0.0;

                forAll(Tc_[zoneI][cellIlocal],i)
                {
                    if(i == 0 || i == (cladSubMeshSize_[zoneI] - 1))
                    {
                        TavClad[cellIglobal] += Tc_[zoneI][cellIlocal][i] * (rcIn_[zoneI] + dC * i) * dC /2.0;
                    }
                    else
                    {
                        TavClad[cellIglobal] += Tc_[zoneI][cellIlocal][i] * (rcIn_[zoneI] + dC * i) * dC;
                    }
                }
                TavClad[cellIglobal] *= 2 / (pow(rcOut_[zoneI],2) - pow(rcIn_[zoneI],2));
            }
        }       
    }
    TavClad.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField>
Foam::constSubscaleFuel::cladOut()
{
    tmp<volScalarField> tToutClad(
        new volScalarField
        (
            IOobject
            (
                "ToutClad",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& ToutClad = tToutClad.ref();
    forAll(Tc_,zoneI)
    {
        forAll(Tc_[zoneI], cellIlocal)
        {
            label zoneId(fuelZoneID_[zoneI]);
            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];  // Mesh cell id for local cellZone id
            ToutClad[cellIglobal] = Tc_[zoneI][cellIlocal][cladSubMeshSize_[zoneI]-1];
        }
    }
    ToutClad.correctBoundaryConditions();
    return tToutClad;
}

// ************************************************************************* //
