scalar TMax(gMax(thermo.T()));
scalar dTMaxdt((TMax-TMaxOld)/runTime.deltaT().value());
Info << "dTMax/dt = " << dTMaxdt << " K/s" << endl;

if ((TMax > TMaxInit) and ((runTime.timeIndex()-runTime.startTimeIndex())>10) and (firstTimePrintTMaxt))
{
    firstTimePrintTMaxt = false;
    Info << "TMax > TMaxInit at t = " << runTime.timeName() << endl;
}

if (TMax >= TBoil)
{
    Info << "Boiling was reached at t = " << runTime.timeName() << endl;
    runTime.writeNow();
    return 0;
}
else if (untilSS)
{
    if ((std::abs(dTMaxdt) < 0.01) and (dTMaxdt > 0))
    {
        SSCounter += 1;
        if (SSCounter >= 5)
        {
            Info << "Steady state was reached after t = " << runTime.timeName() << endl;
            runTime.writeNow();
            return 0;
        }
    }
    else
    {
        SSCounter = 0;
    }            
}

TMaxOld = TMax;

nTimeSteps++;