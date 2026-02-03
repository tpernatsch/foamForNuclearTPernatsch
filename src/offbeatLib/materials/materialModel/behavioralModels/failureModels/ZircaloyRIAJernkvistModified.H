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
    Foam::ZircaloyRIAJernkvistModified

Description
    Overstrain criterion to asses clad failure during the PCMI phase of a RIA transient.

    If the hoop strain is superior to a critical value, the clad is considered failed.
    Critical value of hoop strain is a function of Temperature, strain rate, hydrogen content and fast neutron fluence (E > 1 MeV)

    From: 
    A strain-based clad failure criterion for RIA in light Water Reactors. Lars Olof Jernkvist et al.
    https://inis.iaea.org/collection/NCLCollectionStore/_Public/42/022/42022542.pdf

    The extension from: 
    Investigation of the PCMI failure of pre-hydrided Zircaloy-4 cladding during Reactivity Initiated Accidents 
    with ALCYONE and OFFBEAT fuel performance codes. Matthieu Reymond et al.
    https://www.sciencedirect.com/science/article/pii/S0029549324005302

    is also introduced to take into account the effect of the hydride rim thickness on the hoop plastic strain at failure.
Usage
    In solverDict file:
    \verbatim
    materials
    {
        cladding
        {
            material Zircaloy;

            ...

            failureModel         ZircaloyRiaJernkvistModified;
            meanHydrogenConcentration 400; //mean hydrogen clad concentration
            phi 0.0; //Fast neutron fluence in 1e25 n/m2
            stopIfFailed false;

            rimThicknessReductionFactor true; //- Use additional ductility reduction factor taking into account the rim thickness?
            rimHydrogenConcentration 1200; //rim hydrogen concentration
            rimThickness 80; // rim thickness of the clad, in um
            cladOuterRadius 4.75e-3; // clad outer radius in m
        }

        ...
    }
    \endverbatim

SourceFiles
    ZircaloyRIAJernkvistModified.C
    
\mainauthor
    Matthieu Reymond - CEA/EPFL

\contribution
    J. Sercombe - CEA
    A. Scolaro  - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef ZircaloyRIAJernkvistModified_H
#define ZircaloyRIAJernkvistModified_H

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ZircaloyRIAJernkvistModified Declaration
\*---------------------------------------------------------------------------*/

class ZircaloyRIAJernkvistModified
:
    public failureModel
{
    //- Mean cladding hydrogen concentration (wt.ppm)
    scalar CH_;
    //Fast neutron fluence
    scalar phiValue_;


    //- Use SRecovery ductility recovery factor?
    bool useRimThicknessReductionFactor_;
    //- SReocvery factor parameters
    //- rimThickness (um)
    scalar rimThickness_;
    //- clad outer radius (m)
    scalar RCO_;
    //- rim hydrogen concentration
    scalar CHRim_;


    //- calculated hoopStrainLimit
    volScalarField hoopStrainLimit_;

    //- Cladding Failure Prediction 
    volScalarField CFP_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ZircaloyRIAJernkvistModified(const ZircaloyRIAJernkvistModified&);

        //- Disallow default bitwise assignment
        void operator=(const ZircaloyRIAJernkvistModified&);

protected:
    // protected Member Functions
        virtual scalar computeS0(const scalar&);

public:

    //- Runtime type information
    TypeName("ZircaloyRiaJernkvistModified");


    // Declare run-time constructor selection table
    // Constructors

        //- Construct from dictionary
        ZircaloyRIAJernkvistModified
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~ZircaloyRIAJernkvistModified();


    // Member Functions
    
        //- Return true if material is failed
        virtual bool isFailed(const labelList& addr);

        //- Return name of failure criterion
        virtual word criterionName()
        {
            return typeName_();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
