surfaceScalarField fluxHbyA = fvc::flux(HbyA);

forAll(mesh_.boundary(), patchi)
{
    if (isA<conditionalCyclicAMIFvPatchVectorField>(U.boundaryField()[patchi]))
    {
        const conditionalCyclicAMIFvPatchVectorField& tJunctionPatch =
            refCast< const conditionalCyclicAMIFvPatchVectorField>(U.boundaryField()[patchi]);
        
        fluxHbyA.boundaryFieldRef()[patchi] = fluxHbyA.boundaryFieldRef()[patchi] * tJunctionPatch.fluxCorrectionFactors();

    }
}
