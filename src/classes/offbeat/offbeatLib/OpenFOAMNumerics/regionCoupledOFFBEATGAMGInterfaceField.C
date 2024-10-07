/*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | Website:  https://openfoam.org
     \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
 
 \*---------------------------------------------------------------------------*/
 
 #include "regionCoupledOFFBEATGAMGInterfaceField.H"
 #include "addToRunTimeSelectionTable.H"
 #include "lduMatrix.H"
 
 // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
 namespace Foam
 {
     defineTypeNameAndDebug(regionCoupledOFFBEATGAMGInterfaceField, 0);

     addToRunTimeSelectionTable
     (
         GAMGInterfaceField,
         regionCoupledOFFBEATGAMGInterfaceField,
         lduInterface
     );
     addToRunTimeSelectionTable
     (
         GAMGInterfaceField,
         regionCoupledOFFBEATGAMGInterfaceField,
         lduInterfaceField
     );
#ifdef OPENFOAMESI     
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        regionCoupledOFFBEATGAMGInterfaceField,
        Istream
    );
#endif   
 }
 
 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
 Foam::regionCoupledOFFBEATGAMGInterfaceField::regionCoupledOFFBEATGAMGInterfaceField
 (
     const GAMGInterface& GAMGCp,
     const lduInterfaceField& fineInterface
 )
 :
     GAMGInterfaceField(GAMGCp, fineInterface),
     regionCoupledOFFBEATGAMGInterface_
     (
         refCast<const regionCoupledOFFBEATGAMGInterface>(GAMGCp)
     )
 {}
 
 
 Foam::regionCoupledOFFBEATGAMGInterfaceField::regionCoupledOFFBEATGAMGInterfaceField
 (
    const GAMGInterface& GAMGCp,
#ifdef OPENFOAMESI
    const bool doTransform,
#endif
    const int rank
 )
 :
 #ifdef OPENFOAMFOUNDATION
     GAMGInterfaceField(GAMGCp, rank),
#elif OPENFOAMESI
     GAMGInterfaceField(GAMGCp, doTransform, rank),
#endif     
     regionCoupledOFFBEATGAMGInterface_
     (
         refCast<const regionCoupledOFFBEATGAMGInterface>(GAMGCp)
     )
 {}


#ifdef OPENFOAMESI
Foam::regionCoupledOFFBEATGAMGInterfaceField::regionCoupledOFFBEATGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    Istream& is
)
:
    GAMGInterfaceField(GAMGCp, is),
    regionCoupledOFFBEATGAMGInterface_
    (
        refCast<const regionCoupledOFFBEATGAMGInterface>(GAMGCp)
    )
{}

Foam::regionCoupledOFFBEATGAMGInterfaceField::regionCoupledOFFBEATGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& local,
    const UPtrList<lduInterfaceField>& other
)
:
    GAMGInterfaceField(GAMGCp, local),
    regionCoupledOFFBEATGAMGInterface_
    (
        refCast<const regionCoupledOFFBEATGAMGInterface>(GAMGCp)
    )
{}
#endif
 
 
 // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
 
 Foam::regionCoupledOFFBEATGAMGInterfaceField::~regionCoupledOFFBEATGAMGInterfaceField()
 {}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


#if OPENFOAMESI
void Foam::regionCoupledOFFBEATGAMGInterfaceField::write(Ostream&) const
{}
#endif

// ************************************************************************* //
