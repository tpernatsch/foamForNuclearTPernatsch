/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Info<< "Writing output files\n" << endl;

powerDat << runTime.time().value() << ";" << pTot << endl;

// Write flux optputs

fluxDat << runTime.time().value();

forAll(flux, energyI)
{
    fluxDat << ";" << fvc::domainIntegrate(flux[energyI]).value();
}

fluxDat << endl;

// Write delayed neutron precursors output

precsDat << runTime.time().value();

forAll(prec, precI)
{
    precsDat << ";" << fvc::domainIntegrate(prec[precI]).value();
}

precsDat << endl;

// Write keff

k_effDat << runTime.time().value() << ";" << keff << endl ;

// Write max, min and av fuel and cladding temperature

{

dimensionedScalar Tav = fvc::domainIntegrate(thermo.T()*(TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL))))
                        / max(fvc::domainIntegrate((TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))),dimensionedScalar("",dimensionSet(0,3,0,0,0),SMALL));

dimensionedScalar TavAvFuel = fvc::domainIntegrate(TavFuel)
                        / max(fvc::domainIntegrate((TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))),dimensionedScalar("",dimensionSet(0,3,0,0,0),SMALL));

dimensionedScalar TavAvClad = fvc::domainIntegrate(TavClad)
                        / max(fvc::domainIntegrate((TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))),dimensionedScalar("",dimensionSet(0,3,0,0,0),SMALL));


scalar Tmax = max((thermo.T().internalField())*(TavFuel/(TavFuel+dimensionedScalar("",dimensionSet(0,0,0,1,0),SMALL)))) ;

reduce(Tmax, maxOp<scalar>());

scalar Tmin = min(thermo.T().primitiveField()) ;

reduce(Tmin, minOp<scalar>());

scalar TMaxFuel = max(TinnerFuel.primitiveField()) ;

reduce(TMaxFuel, maxOp<scalar>());

scalar TMaxClad = max(TinnerClad.primitiveField()) ;

reduce(TMaxClad, maxOp<scalar>());

TmaxDat << runTime.time().value() << ";" << Tmax << endl;

TavDat << runTime.time().value() << ";" << Tav.value() << endl;

TavFuelDat << runTime.time().value() << ";" << TavAvFuel.value() << endl;

TavCladDat << runTime.time().value() << ";" << TavAvClad.value() << endl;

TMaxFuelDat << runTime.time().value() << ";" <<  TMaxFuel << endl;

TMaxCladDat << runTime.time().value() << ";" <<  TMaxClad << endl;

TminDat << runTime.time().value() << ";" << Tmin << endl;

totalDisplacedCellsVolumeDat << runTime.time().value() << ";" << (fvc::domainIntegrate((sigmaPow[0]/(sigmaPow[0]+dimensionedScalar("",dimensionSet(1,1,-2,0,0,0,0),SMALL))))).value() << endl;
}

