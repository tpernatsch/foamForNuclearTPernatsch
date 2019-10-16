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

powerDat << runTime.time().value() << ";" << neutronics.power() << endl;

// Write flux optputs

fluxDat << runTime.time().value();

forAll(neutronics.fluxes(), energyI)
{
    fluxDat << ";" << fvc::domainIntegrate(neutronics.fluxes()[energyI]).value();
}

fluxDat << endl;

// Write delayed neutron precursors output

precsDat << runTime.time().value();

forAll(neutronics.precursors(), precI)
{
    precsDat << ";" << fvc::domainIntegrate(neutronics.precursors()[precI]).value();
}

precsDat << endl;

// Write keff

k_effDat << runTime.time().value() << ";" << neutronics.keff() << endl ;

// Write max, min and av fuel and cladding temperature

{
    totalDisplacedCellsVolumeDat 
        << runTime.time().value() << ";" 
        << (
                fvc::domainIntegrate
                (
                        (neutronics.fluxes()[0]
                    /   (
                                neutronics.fluxes()[0]
                            +   dimensionedScalar
                                (
                                    "",
                                    dimensionSet(0,-2,-1,0,0,0,0),
                                    SMALL
                                )
                        )
                    )
                )
            ).value() << endl;
}

