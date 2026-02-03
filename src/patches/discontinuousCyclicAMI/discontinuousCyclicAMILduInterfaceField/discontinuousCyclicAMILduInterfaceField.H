/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::cyclicAMILduInterfaceField

Description
    Abstract base class for cyclic AMI coupled interfaces

SourceFiles
    cyclicAMILduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef discontinuousCyclicAMILduInterfaceField_H
#define discontinuousCyclicAMILduInterfaceField_H

#include "cyclicAMILduInterfaceField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class discontinuousCyclicAMILduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class discontinuousCyclicAMILduInterfaceField
:
    public cyclicAMILduInterfaceField
{
public:

    //- Runtime type information
    TypeName("discontinuousCyclicAMILduInterfaceField");

    //- Constructors
    discontinuousCyclicAMILduInterfaceField() {}

    //- Destructor
    virtual ~discontinuousCyclicAMILduInterfaceField() {}
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
