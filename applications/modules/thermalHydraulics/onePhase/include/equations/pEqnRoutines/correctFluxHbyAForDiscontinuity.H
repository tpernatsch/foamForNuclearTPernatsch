surfaceScalarField fluxHbyA = fvc::flux(HbyA);

forAll(mesh_.boundary(), patchi)
{
    if (isA<scaledDiscontinuousCyclicAMIFvPatchVectorField>(U.boundaryField()[patchi]))
    {
        const scaledDiscontinuousCyclicAMIFvPatchVectorField& tJunctionPatch =
            refCast< const scaledDiscontinuousCyclicAMIFvPatchVectorField>(U.boundaryField()[patchi]);
        
        fluxHbyA.boundaryFieldRef()[patchi] = fluxHbyA.boundaryFieldRef()[patchi] * tJunctionPatch.fluxCorrectionFactors();

    }
}
