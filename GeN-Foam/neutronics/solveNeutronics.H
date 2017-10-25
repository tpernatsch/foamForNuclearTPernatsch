
Info << "Solving neutronics" << endl ;

#include "mapNeutroFields.H"

#include "updateNeutroDiffConst.H"


forAll(flux,energyI)
{
	flux[energyI]=(fluxStar[energyI]-fluxStar2[energyI]*2.0*SP3Neutronics)/discFactor[energyI];
}
if (eigenvalueNeutronics && (runTime.timeIndex()-runTime.startTimeIndex())==1 )
{
	#include "calcPower.H"
	#include "normFluxes.H"
	#include "calcPower.H"
}

#include "initializeNeutroSource.H"
#include "initializeDelayedNeutroSource.H"

if(!eigenvalueNeutronics && couplingIter == 0)
{

	if(integralPredictor)
	{
		#include "integralPredictor.H"
	}
	if(implicitPredictor)	
	{
		#include "implicitPredictor.H"
	}

}	

int aitkenIterCounter = 0;
int aitkenIterNo = 0;
bool aitkenYesNo = true;

int neutroIter = 0;

#include "initializeROMMatrix.H"

#include "createFluxMatrices.H"

do
{
	Info << "Neutron iteration no:  " << neutroIter << nl << endl;

	volFuelPowerNeutro*=0.0;
	neutroResidual = 0.0;

	#include "fluxEq.H"

	if(!solveROM)
	{
		#include "precEq.H"
		#include "initializeDelayedNeutroSource.H"
	}

	neutroIter++;

	Info<< "neutroResidual (couplingIter: " << couplingIter << ",neutroIter: " 
	<< neutroIter << "): " << neutroResidual << " after Aitken no: " << aitkenIterNo << "with solveROM: " << solveROM  << endl; 	

}while(((neutroResidual > neutronIterationResidual) && (neutroIter < maxNeutronIterations) ) &&  !eigenvalueNeutronics);

#include "calcPower.H"

Info<< "Final neutroResidual at coupling/neutroIter iter no " << couplingIter << "/" << neutroIter << ": "  << neutroResidual << " - pow:" << pTot << endl; 

volFuelPowerNeutro /= max(fuelFraction,SMALL) ; //determine vol power IN the fuel

volFuelPower *= 0;
neutroToFluid.mapSrcToTgt(volFuelPowerNeutro, plusEqOp<scalar>(), volFuelPower);
volFuelPower.correctBoundaryConditions();

if (eigenvalueNeutronics)
{
	Info << "pTarget = " << pTarget << nl << endl;
	Info << "keff previous = " << keff << nl << endl;
	keff *= pTot/pTarget;
	keff = min(keff,5);
	keff = max(keff,0.1);
	Info << "keff = " << keff << nl << endl;
	#include "normFluxes.H"
	if (adjustDiscFactors)
	{
		#include "adjustDiscFactors.H"
	}
}

residual = max(residual,neutroResidual);

oneGroupFlux *= 0.0;
forAll(flux,energyI)
{
	oneGroupFlux+=flux[energyI];
}



