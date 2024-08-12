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

Class
    Foam::swellingCorrelationPyC

Description
    Class handling irradiation-induced dimensional change eigenstrain phenomenon
    for PyC with the correlation below for both radial and tangential directions.

    \htmlinclude swelling_swellingCorrelationPyC.html

SourceFiles
    swellingCorrelationPyC.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    January 2023

\*---------------------------------------------------------------------------*/

#ifndef swellingCorrelationPyC_H
#define swellingCorrelationPyC_H

#include "swellingModel.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingCorrelationPyC Declaration
\*---------------------------------------------------------------------------*/

class swellingCorrelationPyC
:
    public swellingModel
{
    // Private data
        //- Name used for the burnup field
        const word fastFluenceName_;

        //- Reference to fast neutron fluence field
        const volScalarField* fastFluence_;

        //- Whether to use spherical coordinate
        bool sp_;
        //- Flux conversion factor. Normally the fast neutron is defined as the neutron
        //- with energy(E) more than 0.1 MeV. However in this model, the neutron (E > 0.18 MeV)
        //- flux is used. So the conversion between flux(>0.1 MeV) to flux(>0.18 MeV) is needed.
        scalar fcf_;

        //- //- Radial and tangential coefficients
        scalarField Ar_;
        scalarField At_;

        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingCorrelationPyC(const swellingCorrelationPyC&);

        //- Disallow default bitwise assignment
        void operator=(const swellingCorrelationPyC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("PyCCorrelation");

    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingCorrelationPyC
        (
            const fvMesh& mesh,
            const dictionary& dict
        );




    // Selectors


    //- Destructor
    virtual ~swellingCorrelationPyC();


    // Member Functions

        //- Update
        virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
