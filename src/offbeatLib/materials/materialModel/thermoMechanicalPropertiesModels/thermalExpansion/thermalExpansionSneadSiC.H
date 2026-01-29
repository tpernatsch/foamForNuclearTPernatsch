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
    Foam::thermalExpansionSneadSiC

Description
    Class modelling the thermal expansion of SiC with a value from Snead. the mean
    thermal expansion coefficient is calculated and the The method of Niffenegger and Reichlin (2012)
    is employed to convert the mean thermal expansion values into instantaneous values.

    Snead model was used for UN TRIOS by Collin(2014).

    -M. Niffenegger and K. Reichlin. The proper use of thermal expansion coefficients
     in finite element calculations. Nuclear Engineering and Design, 243:356–359, 2012.
    -L. L. Snead, T. Nozawa, Y. Katoh, T.-S. Byun, S. Kondo, and D. A. Petti.
     Handbook of sic properties for fuel performance modeling. Journal of Nuclear Materials,
     371:329–377, 2007.
    -Blaise P. Collin. Modeling and analysis of UN TRISO fuel for LWR application
     using the Parfume code. Journal of Nuclear Materials, 451(1):65–77, 2014.


SourceFiles
    thermalExpansionSneadSiC.C

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

#ifndef thermalExpansionSneadSiC_H
#define thermalExpansionSneadSiC_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionSneadSiC Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionSneadSiC
:
    public thermalExpansionModel
{
    // Private data

        //- Snead model parameters
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar Treference;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionSneadSiC(const thermalExpansionSneadSiC&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionSneadSiC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("SiCSnead");

    // Constructors

        //- Construct from dictionary
        thermalExpansionSneadSiC
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel
        );

    //- Destructor
    virtual ~thermalExpansionSneadSiC();


    // Member Functions

    //- Update thermalExpansion
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
