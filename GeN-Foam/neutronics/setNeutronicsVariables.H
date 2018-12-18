
Info<< "Setting neutronics fields\n" << endl;

// Calculate Control Rod position
forAll(CRentries,CRzoneI)
{
    //recalculate positions in case mesh has moved 
    label FOLLOWERzoneId = neutroMesh.cellZones().findZoneID(CRFollowerName[CRzoneI]);

    if(FOLLOWERzoneId == -1)
    {
      Info << "Control rod with follower: " << CRFollowerName(CRzoneI) << " does not exists" << endl;
    }

    scalar cellHeight = -100.0 ;
    scalar maxFollowerHeight = cellHeight ;

    volScalarField centersVSF(neutroMesh.C() & vector(0,0,1));

    Field<scalar> centersSF(centersVSF.internalField());

    forAll(neutroMesh.cellZones()[FOLLOWERzoneId], cellIlocal)
    {

        label cellIglobal = neutroMesh.cellZones()[FOLLOWERzoneId][cellIlocal];

        maxFollowerHeight = max(cellHeight,centersSF[cellIglobal]) ;

        cellHeight = maxFollowerHeight;

    }
    reduce(maxFollowerHeight, maxOp<scalar>());
    CRinitialPosition[CRzoneI] = maxFollowerHeight-initialDistanceFromMeshCR[CRzoneI];
    if(runTime.value()<CRfinish[CRzoneI])
    {
        CRposition[CRzoneI] = (CRinitialPosition[CRzoneI]-CRspeed[CRzoneI]*(max(runTime.value(),CRstart[CRzoneI])-CRstart[CRzoneI]));
    }
    else
    {
        CRposition[CRzoneI] = (CRinitialPosition[CRzoneI]-CRspeed[CRzoneI]*(CRfinish[CRzoneI]-CRstart[CRzoneI]));
    }
    Info << "Control rod " << CRentries[CRzoneI].keyword() << " in position: " << CRposition[CRzoneI] << endl;
}

tmp<volVectorField> tDNeutroRad(new volVectorField(DNeutro));
volVectorField& DNeutroRad = tDNeutroRad.ref();

tmp<volScalarField> tDNeutroAx(new volScalarField(DNeutro & axialOrientation));
volScalarField& DNeutroAx = tDNeutroAx.ref();

DNeutroRad.replace(0,DNeutroRad.component(0)*(1-axialOrientation[0]));
DNeutroRad.replace(1,DNeutroRad.component(1)*(1-axialOrientation[1]));
DNeutroRad.replace(2,DNeutroRad.component(2)*(1-axialOrientation[2]));

DNeutroRad.correctBoundaryConditions();
DNeutroAx.correctBoundaryConditions();

radExp = fvc::div(DNeutroRad)/2.0;
axExp = fvc::grad(DNeutroAx)& axialOrientation;

radExp.correctBoundaryConditions();
axExp.correctBoundaryConditions();

tDNeutroRad.clear();
DNeutroAx.clear();

volScalarField logT( log(max(0.01,(TfuelNeutro/dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef)))) );
volScalarField diffT( (pow(max(TfuelNeutro,dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),SMALL)),0.5)-pow(dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef),0.5)) );

forAll(flux,energyI)
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

        label zoneId = neutroMesh.cellZones().findZoneID(name);

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

            forAll(neutroMesh.cellZones()[zoneId], cellIlocal)
            {
                label cellIglobal = neutroMesh.cellZones()[zoneId][cellIlocal];
                zone = zoneI;

                if(CRNumber>0)
                {
                    if(isFollower && neutroMesh.C()[cellIglobal].z()>CRposition[CRzone])
                    {
                        forAll(entries, zoneII)
                        {
                            if(entries[zoneII].keyword() == CRentries[CRzone].keyword())
                            {
                                zone = zoneII;
                            }
                        }
                    }

                    if(isCR && neutroMesh.C()[cellIglobal].z()<CRposition[CRzone])
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

                D[energyI][cellIglobal] = DList[zone][energyI];
                nuSigmaEff[energyI][cellIglobal] = nuSigmaEffList[zone][energyI];
                sigmaPow[energyI][cellIglobal] = sigmaPowList[zone][energyI];
                sigmaDisapp[energyI][cellIglobal] = sigmaDisappList[zone][energyI];
                forAll(sigmaFromTo,energyJ)
                {
                    sigmaFromTo[energyJ][energyI][cellIglobal] = sigmaFromToList[zone][energyJ][energyI];
                }

                if (parametrize)
                {
                    D[energyI][cellIglobal]
                 += rhoCoolDList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
                  + TCoolDList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempDList[zone][energyI]*logT[cellIglobal]
                      : fuelTempDList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpDList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef)
                  + axialExpDList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpDList[zone][energyI]*(radExp[cellIglobal]);

                    nuSigmaEff[energyI][cellIglobal]
                 += rhoCoolNuSigmaEffList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
                  + TCoolNuSigmaEffList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempNuSigmaEffList[zone][energyI]*logT[cellIglobal]
                      : fuelTempNuSigmaEffList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpNuSigmaEffList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef)
                  + axialExpNuSigmaEffList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpNuSigmaEffList[zone][energyI]*(radExp[cellIglobal]);

                    sigmaPow[energyI][cellIglobal]
                 += rhoCoolSigmaPowList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
                  + TCoolSigmaPowList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempSigmaPowList[zone][energyI]*logT[cellIglobal]
                      : fuelTempSigmaPowList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpSigmaPowList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef)
                  + axialExpSigmaPowList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpSigmaPowList[zone][energyI]*(radExp[cellIglobal]);

                    sigmaDisapp[energyI][cellIglobal]
                 += rhoCoolSigmaDisappList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
                  + TCoolSigmaDisappList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
                  + (
                        fastNeutrons
                      ? fuelTempSigmaDisappList[zone][energyI]*logT[cellIglobal]
                      : fuelTempSigmaDisappList[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpSigmaDisappList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef)
                  + axialExpSigmaDisappList[zone][energyI]*(axExp[cellIglobal])
                  + radialExpSigmaDisappList[zone][energyI]*(radExp[cellIglobal]);

                    forAll(sigmaFromTo,energyJ)
                    {
                        sigmaFromTo[energyJ][energyI][cellIglobal]
                     += rhoCoolSigmaFromToList[zone][energyJ][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
                      + TCoolSigmaFromToList[zone][energyJ][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
                      + (
                            fastNeutrons
                          ? fuelTempSigmaFromToList[zone][energyJ][energyI]*logT[cellIglobal]
                          : fuelTempSigmaFromToList[zone][energyJ][energyI]*diffT[cellIglobal]
                        )
                      + cladExpSigmaFromToList[zone][energyJ][energyI]*(TcladNeutro[cellIglobal] - TcladRef)
                      + axialExpSigmaFromToList[zone][energyJ][energyI]*(axExp[cellIglobal])
                      + radialExpSigmaFromToList[zone][energyJ][energyI]*(radExp[cellIglobal]);
                    }
                }

                if(useGivenDiscFactors || ((runTime.timeIndex()-runTime.startTimeIndex())>0))
                {
                    discFactor[energyI][cellIglobal] = discFactorList[zone][energyI] ;
                }
                else
                {
                    discFactorList[zone][energyI] = discFactor[energyI][cellIglobal] ;
                }

            }//close forall cellIlocal in zoneId
        }//close else from zoneId check
    }//close forall zoneI in zone entries

    D[energyI].correctBoundaryConditions();
    nuSigmaEff[energyI].correctBoundaryConditions();
    sigmaPow[energyI].correctBoundaryConditions();
    sigmaDisapp[energyI].correctBoundaryConditions();
    forAll(sigmaFromTo,energyJ)
    {
        sigmaFromTo[energyJ][energyI].correctBoundaryConditions();
    }
    discFactor[energyI].correctBoundaryConditions();
}

Info << "Done" << endl;
