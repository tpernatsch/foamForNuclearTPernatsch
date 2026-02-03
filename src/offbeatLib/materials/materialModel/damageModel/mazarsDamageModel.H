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
    Foam::mazarsDamageModel

Description
    Damage model for UO2.
    This model describes the loss of stiffness of the material by reducing the elastic modulus
    as the damage increases: E = (1-D)*E.

    The model itself is the so called u-model, based on the work of J. Mazars et al.

    The present model is a parametrization of this model proposed by M. Reymond et al
    to simulate laser heating experiments conducted on UO2. The thermo-mechanical state
    of the ceramic reached during these experiments are close to the ones encountered during a RIA.

    This model is implemented in OFFBEAT for RIA simulations in order to simulate fuel cracking
    and fragmentation during the transient.
    The calculated damage value can also be used to drive the transient fission gas release.

    We thank F. Hamon (EdF R&D) and T. Helfer (CEA) for the numerical implementation available 
    in the MFront software from which the OFFBEAT implementation is largely based on.

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

            rheologyModel mazarsDamageModel;

            rheologyModelOptions
            {

            }
        }

        ...
    }
    \endverbatim

SourceFiles
    mazarsDamageModel.c

\todo 

\mainauthor
    Matthieu Reymond (LRS -- EPFL)

\contribution
    Alessandro Scolaro (LRS -- EPFL)
    Jerome Sercombe (CEA Cadarache)

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef mazarsDamageModel_H
#define mazarsDamageModel_H

#include "damageModel.H"
#include "Table.H"
#include "Time.H"
#include "mechanicsSubSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class mazarsDamageModel Declaration
\*---------------------------------------------------------------------------*/

class mazarsDamageModel
:
    public damageModel
{
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        mazarsDamageModel(const mazarsDamageModel&);

        //- Disallow default bitwise assignment
        void operator=(const mazarsDamageModel&);

protected:

    // Protected data    

        //- Damage
        volScalarField& d_;

        //- Max traction equivalent strain
        volScalarField& Yt_;

        //- Max compression equivalent strain
        volScalarField& Yc_;;
        
        //- Driving variable
        volScalarField& Y_;

        //- Triaxial factor
        volScalarField& r_;

        //- Asymptotic shear
        scalar k_;

        //- Dmax value
        scalar d_max_;

        //- temperatureThreshold activating the crushing damage
        scalar temperatureThreshold_;

        //- Switch for usage of the threshold on equivalent strain rate
        bool equivalentStrainRateThreshold_;

        //- Value of the equivalent srain rate threshold
        scalar strainRateThreshold_;

        //- Do not consider axial cmpt for damage (necessary for smeared
        //  column simulations)
        Switch excludeAxialCmpt_;

        //- Relax factor for damage
        scalar relax_;
        
        //- CurrenTime
        scalar currentTime_;

    // protected Member Functions
        
        //- Compute variable A for traction
        virtual scalar At(const scalar&);

        //- Compute variable A for compression
        virtual scalar Ac(const scalar&);

        //- Compute variable B for traction
        virtual scalar Bt(const scalar&);

        //- Compute variable B for compression
        virtual scalar Bc(const scalar&);

        //- Compute temperature dependent threshold for traction
        virtual scalar epsilonEqt0(const scalar&);

        //- Compute temperature dependent threshold for compression
        virtual scalar epsilonEqc0(const scalar&);

        //- Compute the triaxialFactor for a given strain tensor
        //- Use the solid4FOAM eigen decomposition.
        virtual scalar triaxialFactor(
            const scalar&, const scalar&, const scalar&);

        //- Compute variable A
        virtual scalar A(
            const scalar&, const scalar&, const scalar&, const scalar&);
        
        //- Compute variable B
        virtual scalar B(const scalar&, const scalar&, const scalar&);

public:
 
    //- Runtime type information
    TypeName("mazars");

    // Constructors

        //- Construct from mesh, global options, materials , dict and labelList
        mazarsDamageModel
        (
            const fvMesh& mesh,
            const dictionary& damageDict,
            const dictionary& baseDict
        );


    //- Destructor
    ~mazarsDamageModel();


    // Member Functions

        //- Update damage and apply to stress
        virtual void correct
        (   
            volScalarField& sigmaHyd,
            volSymmTensorField& sigmaDev,
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );

    // Access Functions

        //- Return a constant reference to the damage field
        virtual const volScalarField& d() const
        {
            return d_;
        };
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //