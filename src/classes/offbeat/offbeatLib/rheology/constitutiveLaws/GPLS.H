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
    Foam::GPLS

Description
    Class for the anisotropic viscoplastic model of M. Le Saux for the mechanical behaviour
    of irradiated Zy-4 during a RIA transient.
    
    \htmlinclude constitutiveLaw_GPLS.html

SourceFiles
    GPLS.c

\todo 
Add warning for out of bounds values.
Modify the sigmaEq calculation in the mechanical subSolver to output the Hill's Equivalent Strain and not the Von Mises Equivalent Strain.
It's part of a most larger task consisting in moving stuff like sigma and so on in the rheology class

Remove non necessary outputs 

\mainauthor
    Matthieu Reymond (LRS -- EPFL)

\contribution
    Alessandro Scolaro (LRS -- EPFL)  
    Jerome Sercombe (CEA Cadarache)

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef GPLS_H
#define GPLS_H

#include "constitutiveLaw.H"
#include "Table.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class GPLS Declaration
\*---------------------------------------------------------------------------*/

class GPLS
:
    public constitutiveLaw
{

protected:

    // protected data    
        
        //- fluence value
        scalar phiValue_;

        //- normalised axial position zeta
        scalar zetaValue_;

        //- Viscosity eta scalar field (L*K)
        volScalarField etaVis_;

        //- Relaxtion factor and 
        scalar relax_;

        //theta factor for the semi implicit scheme (not implemented)
        scalar theta_;

        //- tolerance newton raphson
        scalar tol_;

        //- max number of newton rahpson iterations
        scalar maxIter_;

        //- Plastic Strain Tensor
        volSymmTensorField& epsilonP_;        
        
        //- Incremental Plastic Strain tensor
        volSymmTensorField& DEpsilonP_;   
       
        //- Equivalent plastic strain
        volScalarField& epsilonPEq_;

        //- Equivalent plastic strain increment
        volScalarField& DEpsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise 
        volScalarField& activeYield_;

        //- plastic flow direction
        volSymmTensorField& plasticN_;

        //- Maximum allow increase of average epsilon P eq
        scalar maxAverageEpsilonPEq_;

        //- Maximum allowed inscrease of maximum epsilon P Eq
        scalar maxMaximumEpsilonPEq_;

        //- TODO: this variables need more comments
        volScalarField m_;
        volScalarField Kmod_; 
        volScalarField L_;
        volScalarField sigmaEqH_;
        volScalarField tracePlasticN_;
        

    // protected Member Functions
        //- Irradiation damage function
        virtual scalar phiPhi(const scalar&) const;

        //- Strain rate sensivity exponent function and field
        virtual scalar m(const scalar&, const scalar&) const;

        //- Strength coefficient K function and field
        virtual scalar K(const scalar&, const scalar&, const scalar&);

        //- Strain hardening coefficient. Accessories functions + L function + L field
        virtual scalar n0(const scalar&, const scalar&);
        virtual scalar ninf(const scalar&, const scalar&);
        virtual scalar alpha_n(const scalar&);

        virtual scalar L(const scalar&, const scalar&, const scalar&);

        //- Anistropy coefficients functions
        virtual scalar Hrr(const scalar&, const scalar&);
        virtual scalar Hzz(const scalar&, const scalar&);
        virtual scalar Htt(const scalar&, const scalar&);

        //- Hill equivalent strain function + Field
        virtual scalar sigmaEq_func(const symmTensor&, const scalar&, const scalar&);

        //- Contracted double dot product function
        virtual symmTensor normal_tensor(const scalar&, const scalar&, const scalar&, const symmTensor&);

        //- cartesian to cylindrical coordinates
        virtual symmTensor cylSig(const symmTensor&, const vector&);

        //- Function
        virtual scalar F(const scalar&, const scalar&, const scalar&, const scalar&, 
        const scalar&, const scalar&, const scalar&, const scalar&);
        
        //G theta function
        virtual scalar G(const scalar&, const scalar&, const scalar&, const scalar&, const scalar&);

        //- Derivative function for the Newton Raphson loop
        virtual scalar dF(const scalar&, const scalar&, const scalar&, const scalar&, 
        const scalar&, const scalar&, const scalar&, const scalar&);

        //- Disallow default bitwise copy construct
        GPLS(const GPLS&);

        //- Disallow default bitwise assignment
        void operator=(const GPLS&);

public:
 
    //- Runtime type information
    TypeName("GPLS");

    // Constructors

        //- Construct from mesh, global options, materials , dict and labelList
        GPLS
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~GPLS();


    // Member Functions

        //- Integration of the behaviour law and correct fields 
        virtual void correct
        (
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev, 
            volSymmTensorField& epsEl,
            const labelList& addr
        );  

        //- Correct the elastic strain tesnor removing oldTime copmonents of plastic
        // strains  
        virtual void correctEpsilonEl
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );   

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields
        (
            const labelList& addr
        );   

        //- correct the additional strain field
        virtual void correctAdditionalStrain
        (
            volSymmTensorField& additionalStrain,
            const labelList& addr
        );       

        //- Return next time increment according to GPLS model
        virtual scalar nextDeltaT
        (
            const labelList& addr
        ); 

        //- Access functions     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
