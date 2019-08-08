/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "elasticSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
//#include "thermoMechanicsSolver.H"

//#include "gapContactFvPatchVectorField.H"
//#include "gapGasModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elasticSpring::
elasticSpring
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    displacementName_("D")
{}


elasticSpring::
elasticSpring
(
    const elasticSpring& tvtdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tvtdpvf, p, iF, mapper),
    displacementName_(tvtdpvf.displacementName_)
{}


elasticSpring::
elasticSpring
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict),
    displacementName_(dict.lookupOrDefault<word>("displacementName", "D"))
{}


elasticSpring::
elasticSpring
(
    const elasticSpring& tvtdpvf
)
:
    tractionDisplacementFvPatchVectorField(tvtdpvf), 
    displacementName_(tvtdpvf.displacementName_)
{}


elasticSpring::
elasticSpring
(
    const elasticSpring& tvtdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tvtdpvf, iF), 
    displacementName_(tvtdpvf.displacementName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticSpring::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void elasticSpring::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void elasticSpring::updateCoeffs()
{
    if (updated())
    {
        return;
    }
        
    const fvPatchField<vector>& displacement =
        patch().lookupPatchField<volVectorField, vector>("Disp");          
                
    vectorField n(patch().nf());  
    
   /*const fvPatch& patchReference(patch().boundaryMesh().mesh().boundary()[6]);
        
    const scalarField& pressureReference = refCast<const gapContactFvPatchVectorField>
    (
        patchReference.lookupPatchField<volVectorField, vector>("DD")
    ).pressure();        */
   
       const fvPatchField<scalar>& E
        = patch().lookupPatchField<volScalarField, scalar>("E"); 
        
    //const fvPatchField<scalar>& mu =
    //    patch().lookupPatchField<volScalarField, scalar>("muMech");

    //const fvPatchField<scalar>& lambda =
    //    patch().lookupPatchField<volScalarField, scalar>("lambdaMech");        
        
    //word patchName = patch().name();
    //label patchID = patch().boundaryMesh().findPatchID(patchName); 
    //const fvPatchVectorField& totalDispPatch = thermoMechanics.totalDisplacement().boundaryField()[patchID]; 
    
    //- Include the gap gas pressure
   // const gapGasModel& gapGas
   //     = patch().boundaryMesh().mesh().lookupObject<gapGasModel>("gapGas");
      
 
    //traction_ = (E)*(displacement);
    pressure_ = (E)*(displacement&n)/100000;

    //pressure_ = (2*mu+lambda)*(displacement&n)/0.17/10
        
    tractionDisplacementFvPatchVectorField::updateCoeffs();
    

}


void elasticSpring::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    elasticSpring
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

