surfaceScalarField fluxHbyA = fvc::flux(HbyA);
// To be tested if still needed

forAll(mesh_.boundary(), patchi)
{
    if (isA<junctionJumpAMIFvPatchVectorField>(U.boundaryField()[patchi]))
    {
        const junctionJumpAMIFvPatchVectorField& tJunctionPatch =
            refCast< const junctionJumpAMIFvPatchVectorField>(U.boundaryField()[patchi]);
        
        if(tJunctionPatch.discontinuousCyclicAMIPatch().owner())
        {
            fluxHbyA.boundaryFieldRef()[patchi] =tJunctionPatch.scalingFactors()*mesh_.Sf().boundaryField()[patchi] &  (U.boundaryField()[patchi]);
        }
        else
        {
            fluxHbyA.boundaryFieldRef()[patchi] =mesh_.Sf().boundaryField()[patchi] &  (U.boundaryField()[patchi]);      
        }
    }
}
