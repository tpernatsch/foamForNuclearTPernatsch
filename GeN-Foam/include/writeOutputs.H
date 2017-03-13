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

powerDat << runTime.time().value() << "\t" << pTot << endl;

fluxDat << runTime.time().value();

forAll(flux, energyI)
{
    fluxDat << "\t" << fvc::domainIntegrate(flux[energyI]).value();
}

fluxDat << endl;

precsDat << runTime.time().value();

forAll(prec, precI)
{
    precsDat << "\t" << fvc::domainIntegrate(prec[precI]).value();
}

precsDat << endl;

k_eff << keff << endl ;
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

TmaxDat << runTime.time().value() << "\t" << Tmax << endl;

TavDat << runTime.time().value() << "\t" << Tav.value() << endl;

TavFuelDat << runTime.time().value() << "\t" << TavAvFuel.value() << endl;

TavCladDat << runTime.time().value() << "\t" << TavAvClad.value() << endl;

TMaxFuelDat << runTime.time().value() << "\t" <<  TMaxFuel << endl;

TMaxCladDat << runTime.time().value() << "\t" <<  TMaxClad << endl;

TminDat << runTime.time().value() << "\t" << Tmin << endl;

totalDisplacedCellsVolumeDat << runTime.time().value() << "\t" << (fvc::domainIntegrate((sigmaPow[0]/(sigmaPow[0]+dimensionedScalar("",dimensionSet(1,1,-2,0,0,0,0),SMALL))))).value() << endl;
}

