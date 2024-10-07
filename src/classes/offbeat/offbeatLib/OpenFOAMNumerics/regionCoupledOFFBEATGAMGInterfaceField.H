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
    Foam::regionCoupledOFFBEATGAMGInterfaceField
 
Description
    GAMG agglomerated region coupled interface field.
 
SourceFiles
    regionCoupledOFFBEATGAMGInterfaceField.C
 
\*----------------------------------------------------------------------------*/
 
 #ifndef regionCoupledOFFBEATGAMGInterfaceField_H
 #define regionCoupledOFFBEATGAMGInterfaceField_H
 
 #include "GAMGInterfaceField.H"
 #include "regionCoupledOFFBEATGAMGInterface.H"
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                Class regionCoupledOFFBEATGAMGInterfaceField Declaration
 \*---------------------------------------------------------------------------*/
 
 class regionCoupledOFFBEATGAMGInterfaceField
 :
     public GAMGInterfaceField
 {
     // Private Data
 
         //- Local reference cast into the cyclic interface
         const regionCoupledOFFBEATGAMGInterface& regionCoupledOFFBEATGAMGInterface_;
 
 
 public:
 
     //- Runtime type information
     TypeName("regionCoupledOFFBEAT");
 
 
     // Constructors
 
         //- Construct from GAMG interface and fine level interface field
         regionCoupledOFFBEATGAMGInterfaceField
         (
             const GAMGInterface& GAMGCp,
             const lduInterfaceField& fineInterfaceField
         );
 
         //- Construct from GAMG interface and fine level interface field
         regionCoupledOFFBEATGAMGInterfaceField
         (
            const GAMGInterface& GAMGCp,
#ifdef OPENFOAMESI
            const bool doTransform,
#endif
            const int rank
         );


#ifdef OPENFOAMESI
        //- Construct from GAMG interface and fine level interface field
        regionCoupledOFFBEATGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            Istream& is
        );

        //- Construct from GAMG interface and local and remote fields
        regionCoupledOFFBEATGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& local,
            const UPtrList<lduInterfaceField>& other
        );
#endif


         //- Disallow default bitwise copy construction
         regionCoupledOFFBEATGAMGInterfaceField
         (
             const regionCoupledOFFBEATGAMGInterfaceField&
         );

#ifdef OPENFOAMESI
        //- Construct and return a clone
        virtual autoPtr<GAMGInterfaceField> clone
        (
            const GAMGInterface& GAMGCp,
            const UPtrList<lduInterfaceField>& other // other
        ) const
        {
            return autoPtr<GAMGInterfaceField>
            (
                new regionCoupledOFFBEATGAMGInterfaceField(GAMGCp, *this, other)
            );
        }
#endif
 
 
     //- Destructor
     virtual ~regionCoupledOFFBEATGAMGInterfaceField();
 
 
     // Member Functions
 
         // Access
 
             //- Return size
             label size() const
             {
                 return regionCoupledOFFBEATGAMGInterface_.size();
             }
 
 #ifdef OPENFOAMFOUNDATION
             //- Interface matrix update
             virtual void updateInterfaceMatrix
             (
                 scalarField&,
                 const scalarField&,
                 const scalarField&,
                 const direction,
                 const Pstream::commsTypes commsType
             ) const
             {}
#elif OPENFOAMESI            
            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const
            {}


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
#endif             

 
 
     // Member Operators
 
         //- Disallow default bitwise assignment
         void operator=(const regionCoupledOFFBEATGAMGInterfaceField&) = delete;
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //