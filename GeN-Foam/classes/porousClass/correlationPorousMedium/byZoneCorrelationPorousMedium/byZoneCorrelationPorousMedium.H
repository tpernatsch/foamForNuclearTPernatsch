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

Class
    Foam::darcyCorrelationByZonePorousMedium

Description
    Subclass of correlationPorousMedium. Multiple porous media described by correlations.

SourceFiles
    byZoneCorrelationPorousMedium.C

Reference publications
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

#ifndef byZoneCorrelationPorousMedium_H
#define byZoneCorrelationPorousMedium_H

#include "correlationPorousMedium.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class byZoneCorrelationPorousMedium Declaration
\*---------------------------------------------------------------------------*/

class byZoneCorrelationPorousMedium
:
    public correlationPorousMedium
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        byZoneCorrelationPorousMedium(const byZoneCorrelationPorousMedium&);

        //- Disallow default bitwise assignment
        void operator=(const byZoneCorrelationPorousMedium&);


public:

    //- Runtime type information
    TypeName("byZoneCorrelation");


    // Constructors

        //- Construct from mesh
        byZoneCorrelationPorousMedium
        (
            const volScalarField& rho,
            const volVectorField& U,
            const rhoThermo& thermo
        );


    //- Destructor
    virtual ~byZoneCorrelationPorousMedium();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
