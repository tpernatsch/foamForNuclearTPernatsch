forAll(prec,precI)
{
	if(liquidFuel)
	{

		volScalarField diffCoeffNeutro
		(
			IOobject
			(
				"diffCoeffNeutro",
				runTime.timeName(),
				neutroMesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			neutroMesh,
			dimensionedScalar("", dimensionSet(0,2,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		);
		volScalarField diffCoeff = turb.alphat()/thermo.rho()+thermo.mu()/thermo.rho()/ScNo;// (alphaEff=nu/Pr+alphat)
		neutroToFluid.mapTgtToSrc( diffCoeff , plusEqOp<scalar>(), diffCoeffNeutro.primitiveFieldRef());
		diffCoeffNeutro.correctBoundaryConditions();


		fvScalarMatrix precEq
		(
			fvm::ddt(porosityNeutro*(1-eigenvalueNeutronics),prec[precI])
		  	+ fvm::Sp(lambda[precI]*porosityNeutro, prec[precI])
		  	- neutroSource/keff*Beta[precI]
			+ fvm::div(phiNeutro, prec[precI])
			- fvm::laplacian(diffCoeffNeutro,prec[precI])
		        //- fvm::laplacian(dimensionedScalar("", dimensionSet(0,2,-1,0,0,0,0), 0.01),prec[precI])
		        
		);
		precEq.solve();
	}
	else
	{
		fvScalarMatrix precEq
		(
			fvm::ddt((1-eigenvalueNeutronics),prec[precI])
		  	+ fvm::Sp(lambda[precI], prec[precI])
		  	- neutroSource/keff*Beta[precI]
		);
		precEq.solve();
	}
}
