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
    Parent class for slice mapper. This class creates a mapping framework for 
    passing information from an arbitrary OpenFOAM fuel geometry to an idealized 
    1D pin geometry (or viceversa). 
    Essentially the class superimposes a virtual, structured, axial 1D grid 
    on the actual mesh and uses the resulting slices for various operations. 

    Used because:
        - many quantities are calculated on the basis of 1.5D slice average 
          fields (e.g. the relocation needs the slice average burnup or power)
        - when performing 1.5D simulations, it is necessary to perform a slice
          per slice balance of axial forces to correctly take into account the 
          plenum spring and axial friction forces.
        - averaging and min/max operations over slices might be useful for 
          developing output and or plotting toold, mainly for base irradiation
        - might be useful for comparison against standard tools
        - behavioral models that are based on 1-D assumptions and that in the 
        lack of better solutions are still needed in 3-D simulations (e.g.
        classical relocation models)

    The different daugther classes might follow different schemes to split the 
    geometry (e.g. cellSets, by material, etc.). However the basic mechanism 
    remains the same.

    The basic idea is that the slice mapper will provide a mapping between the
    actual 3-D/2-D geometry of OFFBEAT and a division in axial slices,
    superimposed to the geometry. That is, the mapper splits the geometry in
    axial slices and for each slice registers the IDs of the cells that belong
    to that slice.

    The parent class can be selected with the "none" typename.

Usage
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
        const dictionary mapperOptDict_;
        
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
                const dictionary& mapperOptDict
            ),
            (mesh, materials, mapperOptDict)
        );

    // Constructors

        //- Construct from mesh, materials and dict
        sliceMapper
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& mapperOptDict
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
