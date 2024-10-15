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
    Foam::regionCoupledOFFBEATGAMGInterface
    
Description
    GAMG agglomerated coupled region interface.

SourceFiles
    regionCoupledOFFBEATGAMGInterface.C

\*----------------------------------------------------------------------------*/
 
 #ifndef regionCoupledOFFBEATGAMGInterface_H
 #define regionCoupledOFFBEATGAMGInterface_H
 
 #include "regionCoupledBaseOFFBEATGAMGInterface.H"
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 /*---------------------------------------------------------------------------*\
                  Class regionCoupledOFFBEATGAMGInterface Declaration
 \*---------------------------------------------------------------------------*/
 
 class regionCoupledOFFBEATGAMGInterface
 :
     public regionCoupledBaseOFFBEATGAMGInterface
 {
 public:
 
     //- Runtime type information
     TypeName("regionCoupledOFFBEAT");
 
     // Constructors
 
         //- Construct from fine level interface,
         //  local and neighbour restrict addressing
         regionCoupledOFFBEATGAMGInterface
         (
             const label index,
             const lduInterfacePtrsList& coarseInterfaces,
             const lduInterface& fineInterface,
             const labelField& restrictAddressing,
             const labelField& neighbourRestrictAddressing,
             const label fineLevelIndex,
             const label coarseComm
         );
 
         //- Disallow default bitwise copy construction
         regionCoupledOFFBEATGAMGInterface(const regionCoupledOFFBEATGAMGInterface&) = delete;
 
 
     //- Destructor
     virtual ~regionCoupledOFFBEATGAMGInterface();
 
 
     // Member Functions
 
         // I/O
 
             //- Write to stream
             virtual void write(Ostream&) const
             {
                 // TBD. How to serialise the AMI such that we can stream
                 // regionCoupledOFFBEATGAMGInterface.
                 NotImplemented;
             }
 
 
     // Member Operators
 
         //- Disallow default bitwise assignment
         void operator=(const regionCoupledOFFBEATGAMGInterface&) = delete;
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //
