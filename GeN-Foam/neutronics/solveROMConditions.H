

if( ((neutroResidual0-neutroResidual2)/max(SMALL,neutroResidual0)>0.1) )
{
	solveROM = true;
	fullSolutionCounter = 0;
}
else
{
	fullSolutionCounter++;
	if(fullSolutionCounter > 1)
	{
		adjustROM = true;
		ROMAdjustCounter++;
	}
	else
	{
		adjustROM = false;
	}
	solveROM = false;
}

