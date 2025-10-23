surfaceScalarField fluxHbyA = fvc::flux(HbyA);
// To be tested if still needed

forAll(mesh_.boundary(), patchi)
{
    if (isA<branchingJunctionJumpAMIFvPatchVectorField>(U.boundaryField()[patchi]))
    {
        const branchingJunctionJumpAMIFvPatchVectorField& tJunctionPatch =
            refCast< const branchingJunctionJumpAMIFvPatchVectorField>(U.boundaryField()[patchi]);
        
        if(tJunctionPatch.discontinuousCyclicAMIPatch().owner())
        {
            fluxHbyA.boundaryFieldRef()[patchi] =tJunctionPatch.scalingFactors()*mesh_.Sf().boundaryField()[patchi] &  (U.boundaryField()[patchi]);
        }
        else
        {
            fluxHbyA.boundaryFieldRef()[patchi] =mesh_.Sf().boundaryField()[patchi] &  (U.boundaryField()[patchi]);      
        }
    }


    if (isA<connectingJunctionJumpAMIFvPatchVectorField>(U.boundaryField()[patchi]))
    {
        const connectingJunctionJumpAMIFvPatchVectorField& tJunctionPatch =
            refCast< const connectingJunctionJumpAMIFvPatchVectorField>(U.boundaryField()[patchi]);
        const connectingJunctionJumpAMIFvPatchVectorField& nbrtJunctionPatch = 
            refCast< const connectingJunctionJumpAMIFvPatchVectorField>(tJunctionPatch.neighbourPatchField());
        
        if(tJunctionPatch.discontinuousCyclicAMIPatch().owner())
        {
            fluxHbyA.boundaryFieldRef()[patchi] =mesh_.Sf().boundaryField()[patchi] &  ((U.boundaryField()[patchi]));
        }
        else
        {
            fluxHbyA.boundaryFieldRef()[patchi] =nbrtJunctionPatch.scalingFactors()*mesh_.Sf().boundaryField()[patchi] &  ((U.boundaryField()[patchi]));
        }
    }
}
