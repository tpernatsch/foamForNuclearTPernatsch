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
    Foam::sliceMapper

Description
    As the name implies, the `none` axial slice mapper class in OFFBEAT allows you to 
    neglect this aspect of the code, i.e. no virtual slices are defined.

    \note
	The abscence of a mapper might create conflicts with other models that partially
	or entirely depends on it (e.g. some of the behavioral models such as the relocation
	model).
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for selecting the `none`
    mapper class

    In solverDict file:

    \verbatim
        sliceMapper  none;
    \endverbatim

SourceFiles
    sliceMapper.C

\todo
    Remove isFuel?
    Allow different materials inside one slice?

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef sliceMapper_H
#define sliceMapper_H

#include "primitiveFields.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "fvMesh.H"

#include "materials.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sliceMapper Declaration
\*---------------------------------------------------------------------------*/

class sliceMapper
:
    public regIOobject
{
    // Private data
    
protected:

    // Protected data
    
        //- Reference to mesh
        const fvMesh& mesh_;

        // Non const reference to the materials class
        const materials& mat_;

        // Mapper options dictionary
        const dictionary sliceMapperDict_;
        
        //- List of cell IDs for each slice
        mutable PtrList<labelList> sliceAddrList_; 
        
        //- List of cell IDs for each slice (print useful for debugging)
        mutable autoPtr<volScalarField> sliceID_;

        //- List of bool indicating wheter each slice is fuel or not
        mutable List<bool> isFuel_;
        
        //- Number of slices. 
        mutable int nSlices_;
    
    // Protected Member Functions

        // Create sliceID field
        void createSliceID();

        // Update addressing sliceAddrList
        virtual void calcAddressing() const {};

        //- Disallow default bitwise copy construct
        sliceMapper(const sliceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const sliceMapper&);
        
public:

   //- Runtime type information
    TypeName("none");  


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sliceMapper,
            dictionary,
            ( 
                const fvMesh& mesh,
                const materials& materials,
                const dictionary& sliceMapperDict
            ),
            (mesh, materials, sliceMapperDict)
        );

    // Constructors

        //- Construct from mesh, materials and dict
        sliceMapper
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& sliceMapperDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<sliceMapper> New
        (        
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~sliceMapper();


    // Member Functions
    
        //- Return total number of axial slices
        virtual scalar nSlices() const
        {
            return nSlices_;
        } 
    
        //- Return slice addressing    
        virtual const PtrList<labelList>& sliceAddrList() const
        {
            // Topology has changed
            if(mesh_.topoChanging())
            {
                calcAddressing();
            }

            return sliceAddrList_;
        }
        
        //- Return sliceID scalar field (mainly for debugging)
        virtual const volScalarField& sliceID() const
        {
            if(sliceID_.valid())
            {
                return sliceID_();
            }
            else
            {                
                FatalErrorInFunction
                << nl
                << "    Attempting to call sliceID field, but the field has not"
                << " been created."
                << abort(FatalError);

                return volScalarField::null(); // Dummy return field (avoid warning)
            }
        }     
        
        //- Return true/false whether the materials inside the slice are fuel. 
        virtual const List<bool>& isFuel() const
        {
            return isFuel_;
        }   
        
        //- Set and return average of volumetric field passed by argument
        template<class Type>
        const Field<Type>& sliceAverage
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

            
    // IO   
        
        //- Read data from time directory
        virtual bool readData(Istream& is)
        {
                return !is.bad();
        }
        
        //- Write data to time directory
        virtual bool writeData(Ostream& os) const 
        {
                return os.good();
        };    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sliceMapperTemplates.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
