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
    Foam::InterpolateTable

Description
    Class for interpolation in tabulated data.
    This class doesn't actually store any fields, but is intended as
    a quick helper class to perform the actual interpolation for
    different data types, i.e. Fields, FieldFields, etc. effectively 
    facilitating multi-dimensional interpolation.
    Currently provides two interpolation methods, step and linear (see )
    
    This is an alternative to OpenFOAM's built in interpolateTable and
    interpolateTable2D classes, which partially combines the two using
    templates. It was created largely for 2D tables where we want to look up
    one dimension first and then interpolate many times in the resulting
    data, i.e. create a 1D table by interpolating in a 2D table.

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

#ifndef InterpolateTable_H
#define InterpolateTable_H

#include "scalar.H"
#include "scalarField.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Helper base class for interpolation tables
class interpolateTableBase
{
public:
    // Public typedefs
    
        //- Enumeration of interpolation method
        enum interpolationMethod
        {
            STEP,   //!< Each tabulated value applies from that point until 
                    //!  the next point, i.e. step changes.
            LINEAR  //!< Linear interpolation between points
        };

        //- Static lookup of interpolation method names
        static const NamedEnum<interpolationMethod, 2>
            interpolationMethodNames_;
};


/*---------------------------------------------------------------------------*\
                         Class InterpolateTable Declaration
\*---------------------------------------------------------------------------*/

template<class FieldType, class Type, class ReturnType>
class InterpolateTable
:
    interpolateTableBase
{
private:
    
    // Private data

        //- Tabulated points
        const scalarField& xValues_;
        
        //- Tabulated values
        const FieldType& data_;
        
        //- Interpolation method
        interpolationMethod method_;

    // Private Member Functions

        //- Check the table bounds, etc.
        void check() const
        {
            //- Check that table bounds match
            if ((xValues_.size() < 2) || (xValues_.size() != data_.size()))
            {
                //- TODO Extended error message
                FatalErrorInFunction()
                    << "Tabulated data has unexpected/mismatched dimensions:"
                    << endl << tab << "ordinates: " << xValues_.size()
                    << endl << tab << "data:     " << data_.size()
                    << abort(FatalError);
            }
            
            //- Check that supplied x values are in ascending order
            for (int i=1; i < xValues_.size(); i++)
            {
                if (xValues_[i] <= xValues_[i-1])
                {
                    FatalErrorInFunction()
                        << "Supplied table ordinates are not in ascending order."
                        << abort(FatalError);
                }
            }
        }
        
        //- Disallow default bitwise copy construct
        InterpolateTable(const InterpolateTable<FieldType, Type, ReturnType>&);

        //- Disallow default bitwise assignment
        void operator=(const InterpolateTable<FieldType, Type, ReturnType>&);


public:

    // Constructors

        //- Construct from tabulated values
        InterpolateTable
        (
            const scalarField& xValues,
            const FieldType& data,
            const interpolationMethod& method
        )
        :
            xValues_(xValues),
            data_(data),
            method_(method)
        {
            check();
        }


    //- Destructor
    ~InterpolateTable() {}
    
    // Member Operators

        //- Interpolation function
        inline ReturnType operator()(scalar x) const
        {
            if (x >= xValues_[0])
            {
                for (int i=1; i < xValues_.size(); i++)
                {
                    if(x <= xValues_[i])
                    {
                        if (method_ == STEP)
                        {
                            // TODO i-1 or i?
                            return data_[i-1];
                        }
                        else if (method_ == LINEAR)
                        {
                            scalar x0 = xValues_[i-1];
                            scalar x1 = xValues_[i];
                            scalar alpha = (x-x0)/(x1-x0);
                            
                            return (1-alpha)*data_[i-1] + alpha*data_[i];
                        }
                    }
                }
            }
            
            // Out of bounds error
            FatalErrorInFunction()
                << "Supplied value out of bounds." << abort(FatalError);
            
            return 1.0*data_[0];    // Dummy return value
        }

    // Member Functions
    
        //- Return the integral under the tabulated data.
        //- The coefficient k determines the power to which to raise the
        //- x values in calculating the area/volume, e.g. 1 for cartesian (x), 
        //- 2 for cylindrical (r), 3 for spherical (r).
        ReturnType integral(scalar k=1) const
        {
            NotImplemented
            
            return 1.0*data_[0];    // Dummy return value
        }
        
        // Return the next point x_i on the table for the given value of x
        scalar getNextPoint(const scalar x) const
        {
            if (x >= xValues_[0])
            {
                for (int i=1; i < xValues_.size(); i++)
                {
                    if(x < xValues_[i])
                    {
                        return xValues_[i];
                    }
                }
            }
            
            // Out of bounds error
            FatalErrorInFunction()
                << "Supplied value out of bounds." << abort(FatalError);
            
            return GREAT;
        }
        
        // Return the last point x_i on the table
        scalar getLastPoint() const
        {
            return xValues_[xValues_.size()-1];
        }        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
scalar InterpolateTable<scalarField, scalar, scalar>::integral(scalar k) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
