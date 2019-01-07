
Info<< "Setting neutronics fields\n" << endl;

// Calculate Control Rod position
forAll(CRentries,CRzoneI)
{
    //recalculate positions in case mesh has moved 
    label FOLLOWERzoneId = mesh_.cellZones().findZoneID(CRFollowerName[CRzoneI]);

    if(FOLLOWERzoneId == -1)
    {
      Info << "Control rod with follower: " << CRFollowerName(CRzoneI) << " does not exists" << endl;
    }

    scalar cellHeight = -100.0 ;
    scalar maxFollowerHeight = cellHeight ;

    volScalarField centersVSF(mesh_.C() & vector(0,0,1));

    Field<scalar> centersSF(centersVSF.internalField());

    forAll(mesh_.cellZones()[FOLLOWERzoneId], cellIlocal)
    {

        label cellIglobal = mesh_.cellZones()[FOLLOWERzoneId][cellIlocal];

        maxFollowerHeight = max(cellHeight,centersSF[cellIglobal]) ;

        cellHeight = maxFollowerHeight;

    }
    reduce(maxFollowerHeight, maxOp<scalar>());
    CRinitialPosition[CRzoneI] = maxFollowerHeight-initialDistanceFromMeshCR[CRzoneI];
    if(mesh_.time().value()<CRfinish[CRzoneI])
    {
        CRposition[CRzoneI] = (CRinitialPosition[CRzoneI]-CRspeed[CRzoneI]*(max(mesh_.time().value(),CRstart[CRzoneI])-CRstart[CRzoneI]));
    }
    else
    {
        CRposition[CRzoneI] = (CRinitialPosition[CRzoneI]-CRspeed[CRzoneI]*(CRfinish[CRzoneI]-CRstart[CRzoneI]));
    }
    Info << "Control rod " << CRentries[CRzoneI].keyword() << " in position: " << CRposition[CRzoneI] << endl;
}

tmp<volVectorField> tDNeutroRad(new volVectorField(Disp));
volVectorField& DNeutroRad = tDNeutroRad.ref();

tmp<volScalarField> tDNeutroAx(new volScalarField(Disp & axialOrientation));
volScalarField& DNeutroAx = tDNeutroAx.ref();

DNeutroRad.replace(0,DNeutroRad.component(0)*(1-axialOrientation[0]));
DNeutroRad.replace(1,DNeutroRad.component(1)*(1-axialOrientation[1]));
DNeutroRad.replace(2,DNeutroRad.component(2)*(1-axialOrientation[2]));

DNeutroRad.correctBoundaryConditions();
DNeutroAx.correctBoundaryConditions();

volScalarField radExp = fvc::div(DNeutroRad)/2.0;
volScalarField axExp = fvc::grad(DNeutroAx)& axialOrientation;

radExp.correctBoundaryConditions();
axExp.correctBoundaryConditions();

tDNeutroRad.clear();
DNeutroAx.clear();

volScalarField logT( log(max(0.01,(Tfuel/dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef)))) );
volScalarField diffT( (pow(max(Tfuel,dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),SMALL)),0.5)-pow(dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef),0.5)) );

forAll(nuSigmaEff_,energyI)
{

    bool parametrize = true;
    forAll(doNotParametrize, groupI)
    {
        if(energyI == doNotParametrize[groupI])
        {
            parametrize = false;
            Info << "Not Parametrizing XS for group: " << doNotParametrize[groupI] << endl;
        }
    }

    forAll(entries, zoneI)
    {
        label zone = zoneI;

        const word& name = entries[zoneI].keyword();

        Info << "Mesh zone: " << name  << endl;

        label zoneId = mesh_.cellZones().findZoneID(name);

        if(zoneId == -1)
        {
            Info << "WARNING! Mesh zone " << name << " does not exists" << endl;
            Info << "Ignore this warning if it is a completely withdrawn control rod" << endl;
        }
        else
        {
            //check if you are in a follower or CR
            bool isFollower = false;
            bool isCR = false;
            label CRzone = 0;

            forAll(CRentries,CRzoneI)
            {
                const word& CRname = CRentries[CRzoneI].keyword();
                
                if(name == CRname)
                {
                    isCR = true ;
                    CRzone = CRzoneI;
                }
                if(name == CRFollowerName[CRzoneI])
                {
                    isFollower = true ;
                    CRzone = CRzoneI;
                }
            }

            forAll(mesh_.cellZones()[zoneId], cellIlocal)
            {
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
                zone = zoneI;

                if(CRNumber>0)
                {
                    if(isFollower && mesh_.C()[cellIglobal].z()>CRposition[CRzone])
                    {
                        forAll(entries, zoneII)
                        {
                            if(entries[zoneII].keyword() == CRentries[CRzone].keyword())
                            {
                                zone = zoneII;
                            }
                        }
                    }

                    if(isCR && mesh_.C()[cellIglobal].z()<CRposition[CRzone])
                    {
                        forAll(entries, zoneII)
                        {
                            if(entries[zoneII].keyword() == CRFollowerName[CRzone])
                            {
                                zone = zoneII;
                            }
                        }
                    }
                }

                // Adjustment of D and XS for current reactor conditions of:
                // rhoCool, TCool, fuelTemp, cladExp, axialExp and radialExp.

                D_[energyI][cellIglobal] = DList[zone][energyI];
                nuSigmaEff_[energyI][cellIglobal] = nuSigmaEffList[zone][energyI];
                sigmaPow_[energyI][cellIglobal] = sigmaPowList[zone][energyI];
                sigmaDisapp_[energyI][cellIglobal] = sigmaDisappList[zone][energyI];
                forAll(sigmaFromTo_,energyJ)
                {
                    sigmaFromTo_[energyJ][energyI][cellIglobal] = sigmaFromToList[zone][energyJ][energyI];
                }

                if (parametrize)
                {
                    D_[energyI][cellIglobal]
                 += rhoCoolDList[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef)
                  + TCoolDList[zone][energyI]*(TCool[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempDList[zone][energyI]*logT[cellIglobal]
                      : fuelTempDList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpDList[zone][energyI]*(Tclad[cellIglobal] - TcladRef)
                  + axialExpDList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpDList[zone][energyI]*(radExp[cellIglobal]);

                    nuSigmaEff_[energyI][cellIglobal]
                 += rhoCoolNuSigmaEffList[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef)
                  + TCoolNuSigmaEffList[zone][energyI]*(TCool[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempNuSigmaEffList[zone][energyI]*logT[cellIglobal]
                      : fuelTempNuSigmaEffList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpNuSigmaEffList[zone][energyI]*(Tclad[cellIglobal] - TcladRef)
                  + axialExpNuSigmaEffList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpNuSigmaEffList[zone][energyI]*(radExp[cellIglobal]);

                    sigmaPow_[energyI][cellIglobal]
                 += rhoCoolSigmaPowList[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef)
                  + TCoolSigmaPowList[zone][energyI]*(TCool[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempSigmaPowList[zone][energyI]*logT[cellIglobal]
                      : fuelTempSigmaPowList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpSigmaPowList[zone][energyI]*(Tclad[cellIglobal] - TcladRef)
                  + axialExpSigmaPowList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpSigmaPowList[zone][energyI]*(radExp[cellIglobal]);

                    sigmaDisapp_[energyI][cellIglobal]
                 += rhoCoolSigmaDisappList[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef)
                  + TCoolSigmaDisappList[zone][energyI]*(TCool[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempSigmaDisappList[zone][energyI]*logT[cellIglobal]
                      : fuelTempSigmaDisappList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpSigmaDisappList[zone][energyI]*(Tclad[cellIglobal] - TcladRef)
                  + axialExpSigmaDisappList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpSigmaDisappList[zone][energyI]*(radExp[cellIglobal]);

                    forAll(sigmaFromTo_,energyJ)
                    {
                        sigmaFromTo_[energyJ][energyI][cellIglobal]
                     += rhoCoolSigmaFromToList[zone][energyJ][energyI]*(rhoCool[cellIglobal] - rhoCoolRef)
                      + TCoolSigmaFromToList[zone][energyJ][energyI]*(TCool[cellIglobal] - TCoolRef)
                      + (
                            fastNeutrons
                          ? fuelTempSigmaFromToList[zone][energyJ][energyI]*logT[cellIglobal]
                          : fuelTempSigmaFromToList[zone][energyJ][energyI]*diffT[cellIglobal]
                        )
                      + cladExpSigmaFromToList[zone][energyJ][energyI]*(Tclad[cellIglobal] - TcladRef)
                      + axialExpSigmaFromToList[zone][energyJ][energyI]*(axExp[cellIglobal])
                      + radialExpSigmaFromToList[zone][energyJ][energyI]*(radExp[cellIglobal]);
                    }
                }

                if(useGivenDiscFactors || ((mesh_.time().timeIndex()-mesh_.time().startTimeIndex())>0))
                {
                    discFactor_[energyI][cellIglobal] = discFactorList[zone][energyI] ;
                }
                else
                {
                    discFactorList[zone][energyI] = discFactor_[energyI][cellIglobal] ;
                }

            }//close forall cellIlocal in zoneId
        }//close else from zoneId check
    }//close forall zoneI in zone entries

    D_[energyI].correctBoundaryConditions();
    nuSigmaEff_[energyI].correctBoundaryConditions();
    sigmaPow_[energyI].correctBoundaryConditions();
    sigmaDisapp_[energyI].correctBoundaryConditions();
    forAll(sigmaFromTo_,energyJ)
    {
        sigmaFromTo_[energyJ][energyI].correctBoundaryConditions();
    }
    discFactor_[energyI].correctBoundaryConditions();
}

Info << "Done" << endl;
