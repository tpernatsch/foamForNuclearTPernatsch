/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "andersonMixingScheme.H"
#include "scalarMatrices.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar andersonMixingScheme<Type>::projectFields
(
    const Field<Type>& f1, 
    const Field<Type>& f2
) const
{
    scalar v = cmptSum(sumCmptProd(f1, f2));
    reduce(v, sumOp<scalar>());
    return v;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
andersonMixingScheme<Type>::andersonMixingScheme
(
    label size,
    const dictionary& dict
)
:
    accelerationScheme<Type>(size, dict),
    diagonalFactor_(dict.lookupOrDefault<scalar>("diagonalFactor", 1e-4)),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0)),
    order_(readLabel(dict.lookup("order"))),
    x_(order_+1),
    count_(0)
{
    Info<< tab << "alpha: " << alpha_ << nl
        << tab << "order: " << order_ << nl
        << tab << "diagonalFactor: " << diagonalFactor_ << nl
        << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
andersonMixingScheme<Type>::~andersonMixingScheme()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void andersonMixingScheme<Type>::reset()
{
    // Reset everythingl
    count_ = 0;

    forAll(x_, k)
    {
        x_.set(k, NULL);
    }
}


template<class Type>
void andersonMixingScheme<Type>::accelerate()
{
    // Store solution snapshot
    x_.set(count_, new PtrList<volTypeField>(this->psi_.size()));
    
    forAll(this->psi_, k)
    {
        x_[count_].set
        (
            k, 
            new volTypeField(this->psi_[k])
        );
    }
    
    count_++;

    // Apply Anderson mixing if enough snapshots are available
    if (count_ == order_ + 1)
    {
        Info<< "Performing Anderson mixing acceleration" << endl;

        // Build error vectors
        PtrList<FieldField<Field, Type> > e(order_);

        forAll(e, i)
        {
            e.set(i, new FieldField<Field, Type>(this->nFields_));

            forAll(this->psi_, k)
            {
                e[i].set(k, x_[i+1][k].primitiveFieldRef() 
                          - x_[i][k].primitiveFieldRef());
            }
        }

        // Construct projection matrix
        scalarSquareMatrix T(order_, 0.0);
        scalar Tmax = 0.0;

        forAll(e, i)
        {
            for(label j=0; j <= i; j++)
            {
                forAll(this->psi_, k)
                {
                    T[i][j] += projectFields(e[i][k], e[j][k]);
                }

                if (i != j)
                {
                    T[j][i] = T[i][j];
                }
                
                Tmax = max(Tmax, mag(T[i][j]));
            }
        }

        // Renormalise projection matrix
        forAll(e, i)
        {
            for(label j=0; j <= i; j++)
            {
                T[i][j] /= Tmax;
                
                if (i != j)
                {
                    T[j][i] = T[i][j];
                }
            }

            //- Ensure matrix is diagonally dominant
            T[i][i] += diagonalFactor_;
        }
        
        // Solve for optimal coefficients
        scalarField b(order_, 1.0);
        LUsolve(T, b);
        scalar bSum = sum(b);

        if (mag(bSum) < SMALL)
        {
            FatalErrorIn(__FUNCTION__)
                << "Computed coefficients sum up to zero."
                << abort(FatalError);
        }

        b /= bSum;
        
        //- Reconstruct projected solution fields
        forAll(this->psi_, k)
        {
            volTypeField& psi = this->psi_[k];
            psi = dimensioned<Type>("", psi.dimensions(), pTraits<Type>::zero);

            forAll(b, i)
            {
                psi += b[i]*((1-alpha_)*x_[i][k] + alpha_*x_[i+1][k]);
            }
        }

        this->reset();
    }
}


template<class Type>
int andersonMixingScheme<Type>::minIter() const
{
    return order_ + 1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
