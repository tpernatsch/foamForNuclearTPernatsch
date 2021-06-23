/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::regimeMapModel::twoParameters

Description
    Regime map model that depends on two parameters.

Autor
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    twoParameters.C

\*---------------------------------------------------------------------------*/

#ifndef twoParameters_H
#define twoParameters_H

#include "regimeMapModel.H"
#include "regimeDomain2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

/*---------------------------------------------------------------------------*\
                           Class twoParameters Declaration
\*---------------------------------------------------------------------------*/

class twoParameters
:
    public regimeMapModel
{
public:

    enum class interpolationMode
    {
        linear, 
        quadratic
    };

    static const Enum<interpolationMode> interpolationModeNames_;
    
private:

    //- Parameter names
    const word parameter1Name_;
    const word parameter2Name_;

    //- Parameter refs
    const scalarField* parameter1Ptr_;
    const scalarField* parameter2Ptr_;

    //- Table of regime points used to define regimeDomains
    typedef Vector2D<scalar> scalarVector2D;
    typedef HashTable
    <
        scalarVector2D,
        word,
        word::hash
    >   pointTable; 
    pointTable points_;

    //- Table of regimeDomains
    typedef HashTable
    <
        regimeDomain2D,
        word,
        word::hash
    >   regimeDomainTable;
    regimeDomainTable regimeDomains_; 

    //- Scale factors for normalization
    scalar DX_;
    scalar DY_;

    //- Interpolation related
    const interpolationMode interpolationMode_;
    scalar interpolationWidth_;
    bool interpolate_;
    

public:

    //- Runtime type information
    TypeName("twoParameters");

    //- Constructors

        twoParameters
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~twoParameters(){}

    //- Member Functions

        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //