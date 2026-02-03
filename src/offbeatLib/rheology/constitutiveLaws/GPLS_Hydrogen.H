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
    Foam::hydrogenGPLS

Description
    hydrogenGPLS (Guisson-Poussard-Le Saux) mechanical behavioral law
    Viscoplastic behavior for fresh and irradiated zyr4 under RIA loading conditions. Used in ALCYONE (CEA - FRANCE) fuel code for RIA simulations. NOT for base irradiation.
    Elastic behavior is assumed isotropic and is given by Hooke's Law.
    Viscoplastic behavior takes into account the orthotropic nature of irradiated ZYR. Make uses of the Hill Tensor and Hill's equivalent strain.
    No stress threshold between elastic and viscoplastic flow.

    The fluence (phi) is required as an input in 1e25 n.m-2 to take into account irradiation effect. IE., put 1 for 1e25

    A second (non physical) parameter zeta is needed.
    In the model, the parameter zeta corresponds to the normalised axial position in the (total) fuel rod (0= bottom, 1=top)
    Introduced to take into account the temperature axial profile during base irradiation
    RIA test pulses in research reactor are typically done on refabricated fuels elements previously irradiated in commercial reactors 
    The strength modulus K is weakly dependant on this parameter
    Default value is 0.5

    Validity range:
    fluence/BU up to 10e25 .nm-2 / 64 GWd/t(U)
    20°C < T < 1100 °C ---- (in practice, in RIA simulations, we can go over 1100 °C)
    3e-4 < ddt(epsilon) < 5 .s-1
    
    from:
    A model to describe the anisotropic viscoplastic mechanical behavior of fresh and irradiated Zircaloy-4 fuel claddings under RIA loading conditions
    Journal of Nuclear Materials, vol 378.1, 60-69, 2008
    https://doi.org/10.1016/j.jnucmat.2008.04.017
    
    Le Saux thesis (french - english) : https://www.theses.fr/2008ENMP1608
    
    Interested reader may also refers to:
    The PROMETRA Program: Fuel Cladding Mechanical Behavior under High Strain Rate
    Nuclear Technology, vol 157.3, 215-229, 2007
    https://doi.org/10.13182/NT07-A3814

Usage
    \verbatim
    rheology byMaterial;
    
    ...

    // List of materials, one per cellZone.
    materials
    {
        ...

        cladding
        {
            material Zircaloy;
            ...

            rheologyModel hydrogenGPLS;

            rheologyModelOptions
            {
                phi      1;
                zeta    0.5;

                error 1e-6; //- Optionnal. Tolerance for exiting the newton-rahpson iterations (default is 1e-6)
                max_iter 10000; //- Optional. Max number of newton-rahpson iterations (default is 10 000)
                relax 0.01; //- Optional. Relaxation factor for newton raphson iterations (default is 0.01)
                theta 1; //- Optional. Theta value (between 0 and 1) for the integration scheme. 0 --> fully explicit (Forward Euler). 
                                                                                                 1 --> fully implicit (Backard Euler)
                                                                                                 0.5 --> semi-implicit (mid point rule / Cranck-Nicholson)
                Fully explicit is not advised. Use 1 or 0.5. The latter is the most efficient and robust scheme.
                More infos at https://tfel.sourceforge.net/documents/mfront/behaviours.pdf (french).
            }
        }

        ...
    }
    \endverbatim

SourceFiles
    hydrogenGPLS.c

\todo Clean variables and parameters as some are not useful anymore
      Add warning for out of bounds values.
      Modify the sigmaEq calculation in the mechanical subSolver to output the Hill's Equivalent Strain and not the Von Mises Equivalent Strain.
      It's part of a most larger task consisting in moving stuff like sigma and so on in the rheology class
      Tangent matrix operator? Solve on both dp and dEpsEl? Check if it is sensible to do so. (Jacobian and jacobian inversion, mainly)

\mainauthor
    Matthieu Reymond (LRS -- EPFL)

\contribution

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef hydrogenGPLS_H
#define hydrogenGPLS_H

#include "constitutiveLaw.H"
#include "Table.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hydrogenGPLS Declaration
\*---------------------------------------------------------------------------*/

class hydrogenGPLS
:
    public constitutiveLaw
{

protected:

    // protected data    
        
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
        
        //- Concentration of H in solid solution
        volScalarField& Css_;

        //- Concentration of H in solid hydrides
        volScalarField& Cpp_;

        //- Concentration of H total
        volScalarField& Ctot_;

        //- Strain rate sensitivity exponenent
        volScalarField m_;

        //- Strength coefficient
        volScalarField Kmod_; 

        //- Strain hardening coefficient
        volScalarField L_;

        //- Viscosity eta scalar field (L*K)
        volScalarField etaVis_;

        //- Equivalent Hill stress
        volScalarField sigmaEqH_;

        //- Maximum allow increase of average epsilon P eq
        scalar maxAverageEpsilonPEq_;

        //- Maximum allowed inscrease of maximum epsilon P Eq
        scalar maxMaximumEpsilonPEq_;
        
    // protected Member Functions
        //- Strain rate sensivity exponent function
        virtual scalar m(const scalar&);

        //- Strength coefficient K function
        virtual scalar K(const scalar&, const scalar&, const scalar&);

        //- Strain hardening coefficient. Accessories functions + L function + L field
        virtual scalar n0(const scalar&, const scalar&);
        virtual scalar alpha_n(const scalar&, const scalar&);

        virtual scalar L(const scalar&, const scalar&, const scalar&);

        //- Anistropy coefficients functions
        virtual scalar Hrr(const scalar&);
        virtual scalar Hzz(const scalar&);
        virtual scalar Htt(const scalar&);

        //- Hill equivalent strain function + Field
        virtual scalar sigmaEq_func(const symmTensor&, const scalar&);

        //- Contracted double dot product function
        virtual symmTensor normal_tensor(const scalar&, const scalar&, const symmTensor&);

        //- cartesian to cylindrical coordinates
        virtual symmTensor cylSig(const symmTensor&, const vector&);

        //- Function
        virtual scalar F(const scalar&, const scalar&, const scalar&, 
                         const scalar&, const scalar&, const scalar&, 
                         const scalar&, const scalar&);
        
        //G theta function
        virtual scalar G(const scalar&, const scalar&, const scalar&, 
                         const scalar&, const scalar&);

        //- Derivative function for the Newton Raphson loop
        virtual scalar dF(const scalar&, const scalar&, const scalar&,
						  const scalar&, const scalar&, const scalar&,
                          const scalar&, const scalar&);

        void correctHydrogenConcentrations(const labelList&);

        //- Disallow default bitwise copy construct
        hydrogenGPLS(const hydrogenGPLS&);

        //- Disallow default bitwise assignment
        void operator=(const hydrogenGPLS&);

public:
 
    //- Runtime type information
    TypeName("hydrogenGPLS");

    // Constructors

        //- Construct from mesh, global options, materials , dict and labelList
        hydrogenGPLS
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~hydrogenGPLS();


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

        //- Return next time increment according to hydrogenGPLS model
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
