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
    Foam::swellingPARFUMEBuffer

Description
    Class handling irradiation-induced dimensional change phenomenon for Buffer
    derived from PARFUME manual.

    \htmlinclude swelling_swellingPARFUMEBuffer.html

SourceFiles
    swellingPARFUMEBuffer.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    December 2022

\*---------------------------------------------------------------------------*/

#ifndef swellingPARFUMEBuffer_H
#define swellingPARFUMEBuffer_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"
#include "Function1.H"
#include "Table.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingPARFUMEBuffer Declaration
\*---------------------------------------------------------------------------*/

class swellingPARFUMEBuffer
:
    public swellingModel
{
    // Private data
        //- Name used for the density field
        const word densityName_;

        //- Reference to density field
        const volScalarField* rho_;

        //- Name used for the fast neutron fluence field
        const word fastFluenceName_;

        //- Reference to fast neutron fluence field
        const volScalarField* fastFluence_;

        //- Name used for the fast neutron flux field
        const word fastFluxName_;

        //- Reference to fast neutron flux field
        const volScalarField* fastFlux_;
        //- Flux conversion factor. Normally the fast neutron is defined as the neutron
        //- with energy(E) more than 0.1 MeV. However in this model, the neutron (E > 0.18 MeV)
        //- flux is used. So the conversion between flux(>0.1 MeV) to flux(>0.18 MeV) is needed.
        scalar fcf_;

        //- Perturbation parameter
        scalar perturb;

        interpolateTableBase::interpolationMethod Method_;
        interpolateTableBase::outofBoundsMehtod obMethod_ex_;
        interpolateTableBase::outofBoundsMehtod obMethod_fix_;

        scalarField tempValues_;
        scalarField BAFValues_r_;
        scalarField BAFValues_t_;
        scalarField densityValue_;
        scalarField strain_iso_;

        FieldField<Field, scalar> a1_r_;
        FieldField<Field, scalar> a2_r_;
        FieldField<Field, scalar> a3_r_;
        FieldField<Field, scalar> a4_r_;

        autoPtr<scalarFieldInterpolateTable> table_a1_r_;
        autoPtr<scalarFieldInterpolateTable> table_a2_r_;
        autoPtr<scalarFieldInterpolateTable> table_a3_r_;
        autoPtr<scalarFieldInterpolateTable> table_a4_r_;

        autoPtr<scalarInterpolateTable> table_epsilon_iso_;




    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingPARFUMEBuffer(const swellingPARFUMEBuffer&);

        //- Disallow default bitwise assignment
        void operator=(const swellingPARFUMEBuffer&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("BufferPARFUME");

    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingPARFUMEBuffer
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors


    //- Destructor
    virtual ~swellingPARFUMEBuffer();


    // Member Functions

        //- Update
        virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
