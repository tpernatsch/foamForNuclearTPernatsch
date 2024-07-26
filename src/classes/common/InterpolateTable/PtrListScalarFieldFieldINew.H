/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::PtrListScalarFieldFieldINew

Description
    Constructor class for PtrList<FieldField<Field, scalar>>

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE 
    DE LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems 
    Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE 
    DE LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems 
    Behaviour)


\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef PtrListScalarFieldFieldINew_H
#define PtrListScalarFieldFieldINew_H

#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
//- Constructor class for FieldField<Field, scalar>
class PtrListScalarFieldFieldINew
{
    public:
        
        //- Construct from istream
        PtrListScalarFieldFieldINew() {}
        
        //- New from input stream
        autoPtr<FieldField<Field, scalar>> operator()(Istream& is) const
        {
            return autoPtr<FieldField<Field, scalar>>(new FieldField<Field, scalar>(PtrList<scalarField>(is, scalarFieldFieldINew())));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
