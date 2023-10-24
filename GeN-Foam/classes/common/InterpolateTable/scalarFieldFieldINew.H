/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::scalarFieldFieldINew

Description
    Constructor class for FieldField<Field, scalar>

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

#ifndef scalarFieldFieldINew_H
#define scalarFieldFieldINew_H

#include "scalarFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
//- Constructor class for FieldField<Field, scalar>
class scalarFieldFieldINew
{
    public:
        
        //- Construct from istream
        scalarFieldFieldINew() {}
        
        //- New from input stream
        autoPtr<scalarField> operator()(Istream& is) const
        {
            return autoPtr<scalarField>(new scalarField(is));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
