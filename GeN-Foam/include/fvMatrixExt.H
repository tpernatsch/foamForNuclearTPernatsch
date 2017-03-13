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

Class
    Foam::subscaleFuel

Description
    Base class used to calculate behavior of fuel. Details of fuel are specified in subclasses.
    It uses a dictionary called fuelProperties

Reference publications
	Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
	"GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
	analysis of nuclear reactors", Nuclear Engineering and Design, submitted

	Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
	solver to the European Sodium Fast Reactor and verification against available codes",
	Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226

Extended fvMatrix to have matrix multiplication and source with effect of boundaries and including already dividion by cell volume. Works only for scalars!

Author
    Carlo Fiorina <carlo.fiorina@hotmail.it; carlo.fiorina@psi.ch; carlo.fiorina@polimi.it; 
                   carlo.fiorina@epfl.ch;>


\*---------------------------------------------------------------------------*/

#ifndef fvMatrixExt_H
#define fvMatrixExt_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class fvMatrixExt Declaration
\*---------------------------------------------------------------------------*/
template<class Type>
class fvMatrixExt
:
    public fvMatrix<Type>
{

public:


// Constructors

// Construct from mesh and dictionary
fvMatrixExt(const tmp<fvMatrix<Type> >& tfvm):
	fvMatrix<Type>(tfvm)
{}

//- Destructor
~fvMatrixExt()
{}

// Member Functions


Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > Source() const
{


    tmp<GeometricField<Type, fvPatchField, volMesh> > tSource
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "Source("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Source = tSource();

    Source.primitiveFieldRef() = this->source();
    Source.primitiveFieldRef() /= this->psi().mesh().V();
    Source.correctBoundaryConditions();

    return tSource;
		
}


Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > SourceWithBC() const //beta version, may give problems
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tSource
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "Source("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Source = tSource.ref();

    Source.primitiveFieldRef() = this->source();
    Foam::fvMatrix<Type>::addBoundarySource(Source.primitiveFieldRef(),false);//false exclude processor boundaries
    Source.primitiveFieldRef() /= this->psi().mesh().V();
    Source.correctBoundaryConditions();	
    return tSource;
		
}


Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBC() 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), this->psi().internalField(), this->boundaryCoeffs(),  this->psi().boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(this->psi().size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*this->psi().internalField();
    AmulWithBC.primitiveFieldRef() /= this->psi().mesh().V();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;		
}



Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBC(GeometricField<Type, fvPatchField, volMesh>& phi) 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+phi.name()+')',
                phi.instance(),
                phi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi.mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), phi.internalField(), this->boundaryCoeffs(),  phi.boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(phi.size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*phi.internalField();
    AmulWithBC.primitiveFieldRef() /= phi.mesh().V();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;		
}




Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > SourceWithBCnoV() const //beta version, may give problems
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tSource
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "Source("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions(),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Source = tSource.ref();

    Source.primitiveFieldRef() = this->source();
    Foam::fvMatrix<Type>::addBoundarySource(Source.primitiveFieldRef(),false);//false exclude processor boundaries
    Source.correctBoundaryConditions();
    return tSource;
		
}


Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBCnoV() 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions(),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), this->psi().internalField(), this->boundaryCoeffs(),  this->psi().boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(this->psi().size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*this->psi().internalField();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;
	
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBCnoV(GeometricField<Type, fvPatchField, volMesh>& phi) 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+phi.name()+')',
                phi.instance(),
                phi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi.mesh(),
            this->dimensions(),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), phi.internalField(), this->boundaryCoeffs(),  phi.boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(phi.size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*phi.internalField();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;			
}



Type initialResidual() 
{
    scalarField saveDiag(this->diag());

    Field<Type> Source(this->source());

    Foam::fvMatrix<Type>::addBoundarySource(Source,true);//false exclude processor boundaries

     scalarField psi = this->psi().primitiveField();
    Foam::fvMatrix<Type>::addBoundaryDiag(this->diag(), 0);

         FieldField<Field, scalar> bouCoeffs
         (
             this->boundaryCoeffs()
         );
 
         FieldField<Field, scalar> intCoeffs
         (
             this->internalCoeffs()
         );
 
         lduInterfaceFieldPtrsList interfaces =
             this->psi().boundaryField().scalarInterfaces();


         this->initMatrixInterfaces
         (
             bouCoeffs,
             interfaces,
             psi,
             Source,
             0
         );
 
         this->updateMatrixInterfaces
         (
             bouCoeffs,
             interfaces,
             psi,
             Source,
             0
         );


     label nCells = psi.size();

     scalarField wA(nCells);
     scalarField pA(nCells);

    this->Foam::fvMatrix<Type>::Amul(wA, psi, this->boundaryCoeffs(),  this->psi().boundaryField().scalarInterfaces() ,  direction(0));


     // --- Calculate initial residual field
     scalarField rA(Source - wA);

     // --- Calculate normalisation factor

this->Foam::fvMatrix<Type>::sumA(pA, bouCoeffs, interfaces);

pA *= gAverage(psi,this->Foam::lduMatrix::mesh().comm());

//this->Foam::lduMatrix::mesh().comm()

scalar normFactor = 
gSum
(
       (mag(wA - pA) + mag(Source - pA))(),
       this->Foam::lduMatrix::mesh().comm()
)
+ SMALL;


    this->diag() = saveDiag;

    return gSumMag(rA,this->Foam::lduMatrix::mesh().comm())  / normFactor;

		
}



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
