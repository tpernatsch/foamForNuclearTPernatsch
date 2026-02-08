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
    Foam::swellingParfumePyC

Description
    Class handling irradiation-induced dimensional change eigenstrain phenomenon
    for PyC derived from Parfume manual. This class is the general situation for
    the class: swellingParfumeBuffer (BAF=1).

    \htmlinclude swelling_swellingParfumePyC.html

SourceFiles
    swellingParfumePyC.C

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

#ifndef swellingParfumePyC_H
#define swellingParfumePyC_H

#include "swellingModel.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingParfumePyC Declaration
\*---------------------------------------------------------------------------*/

class swellingParfumePyC
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

        //- The as-fabricated anisotropy (dimensionless)
        scalar BAF_;
        //- Whether to use spherical coordinate
        bool sp_;
        //- Flux conversion factor. Normally the fast neutron is defined as the neutron
        //- with energy(E) more than 0.1 MeV. However in this model, the neutron (E > 0.18 MeV)
        //- flux is used. So the conversion between flux(>0.1 MeV) to flux(>0.18 MeV) is needed.
        const scalar fcf_;

        interpolateTableBase::interpolationMethod Method_;
        interpolateTableBase::outOfBoundsMethod obMethod_ex_;
        interpolateTableBase::outOfBoundsMethod obMethod_fix_;

        scalarField tempValues_;
        scalarField BAFValues_r_;
        scalarField BAFValues_t_;
        scalarField densityValue_;
        scalarField strain_iso_;
        scalarField strain_r_t_;

        FieldField<Field, scalar> a1_r_;
        FieldField<Field, scalar> a2_r_;
        FieldField<Field, scalar> a3_r_;
        FieldField<Field, scalar> a4_r_;
        FieldField<Field, scalar> a1_t_;
        FieldField<Field, scalar> a2_t_;
        FieldField<Field, scalar> a3_t_;
        FieldField<Field, scalar> a4_t_;

        autoPtr<scalarFieldInterpolateTable> table_a1_r_;
        autoPtr<scalarFieldInterpolateTable> table_a2_r_;
        autoPtr<scalarFieldInterpolateTable> table_a3_r_;
        autoPtr<scalarFieldInterpolateTable> table_a4_r_;
        autoPtr<scalarFieldInterpolateTable> table_a1_t_;
        autoPtr<scalarFieldInterpolateTable> table_a2_t_;
        autoPtr<scalarFieldInterpolateTable> table_a3_t_;
        autoPtr<scalarFieldInterpolateTable> table_a4_t_;

        autoPtr<scalarInterpolateTable> table_epsilon_iso_;
        autoPtr<scalarInterpolateTable> table_epsilon_r_t_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingParfumePyC(const swellingParfumePyC&);

        //- Disallow default bitwise assignment
        void operator=(const swellingParfumePyC&);

    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(swellingParfumePyC, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(swellingParfumePyC, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(swellingParfumePyC, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("PyCParfume");

    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingParfumePyC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );




    // Selectors


    //- Destructor
    virtual ~swellingParfumePyC();


    // Member Functions

        //- Update
        virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
