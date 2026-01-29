/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gapTRISO.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "fuelMaterial.H"
#include "unitConversion.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gapTRISO, 0);
    addToRunTimeSelectionTable
    (
        gapGasModel,
        gapTRISO,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::gapTRISO::calcInitialMass()
{
    //- Total volume available to the gas
    const scalar Vtot =
    gapV_ + gapOffsetV_ + bufferV_ + (reserveV_) + crackV_;

    //- Volume weighted gas temperature in the rod
    const scalar gasT =
    (
        gapV_*gapT_ + gapOffsetV_*gapT_ + bufferV_*bufferT_
        + reserveV_*reserveT_ + crackV_*crackT_)/Vtot;

    //- Gas density using ideal gas law
    scalar rho = rhoMixture(gasP_, gasT);

    //- Gas initial mass, in kg
    gasM0_ = rho*Vtot ;

    //- Update gas mass
    gasM_ = gasM0_;
}


Foam::scalar Foam::gapTRISO::rhoMixture
(
    const scalar& gasP,
    const scalar& gasT
) const
{
    scalar rhoMixture(0.0);
    scalar R_mixture(0.0);
    scalar R = Foam::constant::physicoChemical::R.value();

    forAll(species_, i)
    {
        R_mixture += R/speciesW_[i]*1000*Y_[i];
    }

    if (R_mixture == 0.0)
    {
      rhoMixture = 0.0;
    }
    else
    {
      rhoMixture = gasP/(R_mixture*gasT);
    }

    return rhoMixture;
}


void Foam::gapTRISO::findPatchIDs(const word groupName)
{
    if(!patchIDsPtr_.valid())
    {
        patchIDsPtr_.reset(new HashTable<std::set<label>, word>(10));
    }

    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();

    //- Find associated patches on the mesh that bound the gap gas
    if
    (
#ifdef OPENFOAMFOUNDATION        
        !gapGasDict_.isNull() 
#elif OPENFOAMESI
        !gapGasDict_.isNullDict() 
#endif        
        and gapGasDict_.found(groupName+"Patches")
    )
    {
        wordList gapPatches(gapGasDict_.lookup(groupName+"Patches"));

        forAll(gapPatches, patchI)
        {
            word gapPatchName(gapPatches[patchI]);

            label gapPatchID
            (mesh_.boundaryMesh().findPatchID(gapPatchName));

            if(gapPatchID == -1)
            {
                FatalErrorIn("Foam::gapTRISO::findPatchIDs(const word groupName)")
                << gapPatchName << " patch from list " << groupName << "Patches"
                << " not found in mesh." << exit(FatalError);
            }

            HashTable<std::set<label>, word>::iterator iter = patchIDs.find
            (
                groupName
            );

            if (iter != patchIDs.end())
            {
                iter().insert(gapPatchID);
            }
            else
            {
                std::set<label> s;
                s.insert(gapPatchID);
                patchIDs.insert(groupName, s);
            }
        }

        if(!gapPatches.size())
        {
            //- Insert an empty set --> allow an empty patch list
            std::set<label> s;
            patchIDs.insert(groupName, s);
        }
    }

    const Foam::HashTable<Foam::labelList, Foam::word>& groupPatchIDs =
    mesh_.boundaryMesh().groupPatchIDs();

    if(groupPatchIDs.find(groupName) != groupPatchIDs.end())
    {
        labelList gapPatchIDs(mesh_.boundaryMesh().groupPatchIDs()[groupName]);

        forAll(gapPatchIDs, patchI)
        {
            label gapPatchID(gapPatchIDs[patchI]);

            HashTable<std::set<label>, word>::iterator iter = patchIDs.find
            (
                groupName
            );

            if (iter != patchIDs.end())
            {
                iter().insert(gapPatchID);
            }
            else
            {
                std::set<label> s;
                s.insert(gapPatchID);
                patchIDs.insert(groupName, s);
            }
        }
    }

    if
    (
        patchIDs.find(groupName) == patchIDs.end()
    )
    {
        FatalErrorIn("Foam::gapTRISO::findPatchIDs(const word groupName)")
        << groupName << " patches not found." << nl
        << "Either: " << nl
        << " - group the patches using the group name \""
        << groupName << "\" or" << nl
        << " - list the patch names in \"" << groupName << "Patches\" "
        << "in the gapGasOptions dict in solverDict file"  << nl
        << " (the list can be empty)"  << nl ;
        Info << exit(FatalError);
    }
}

void Foam::gapTRISO::correctGap()
{
    //- Initialize useful quantities
    scalar volume(0);
    scalar area(0);
    scalar temperature(0.0);
    scalar VoverT(0.0);
    //- Loop over specified patches to calculate total gas volume from the
    //- geometry. If the cladding surface is larger than the fuel surface, the
    //- calculated volume will include the plena.
    HashTable<std::set<label>, word>& patchIDs = patchIDsPtr_();
    std::set<label>& gapPatchIDs(patchIDs["gap"]);

    for(const label& patchID : gapPatchIDs)
    {
        const fvPatch& p = mesh_.boundary()[patchID];

        const regionCoupledOFFBEATFvPatch& patch
            = refCast<const regionCoupledOFFBEATFvPatch>(p);

        const regionCoupledOFFBEATFvPatch& nbrPatch
            = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

        //- Current patch quantities

        //- Patch center face vectors, displacement and temperature fields
        const vectorField& cf = patch.Cf();
        const fvPatchVectorField& Df = DorDD_.boundaryField()[patchID];
        const fvPatchScalarField& T = T_.boundaryField()[patchID];

        //- Patch normal vectors
        // TODO: apply F to rotate normal in case of large strain
        const vectorField nf = patch.nf();


        //- Take face areas from AMI as it is more precise than the patch
        // (only if updateAMI is true, otherwise it gives the same value)
        const AMIInterpolation& ami = patch.owner() ?
        patch.regionCoupledPatch().AMI()
        :
        patch.nbrPatch().regionCoupledPatch().AMI();

        List<scalar> magSf = patch.owner() ? 
        ami.srcMagSf() : ami.tgtMagSf();

        //- Scale the area with surface overlap scaling factors
        scalarField factors = scalingFactors_[patchID];
        magSf = magSf*factors;

        //- Neihgbor patch quantities

        const word nbrPatchName = nbrPatch.name();
        const label nbrPatchID =
        mesh_.boundaryMesh().findPatchID(nbrPatchName);

        const vectorField nbrCf =
        patch.regionCoupledPatch().interpolate(nbrPatch.Cf());

        const vectorField nbrDf =
        patch.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volVectorField, vector>(DorDD_.name())
        );

        const vectorField nbrNf =
        patch.regionCoupledPatch().interpolate(nbrPatch.nf());

        //- The nbr face areas are obtained by scaling the owner areas by the
        //  fraction of target face covered by source faces
        scalarField nbrScaleFactor = patch.owner() ?
        patch.regionCoupledPatch().interpolate
        (
            ami.tgtWeightsSum()/max(scalingFactors_[nbrPatchID], VSMALL)
        )
        :
        scalingFactors_[patchID]/max(ami.tgtWeightsSum(), VSMALL);

        scalarField nbrMagSf = magSf/max(nbrScaleFactor, VSMALL);

        //- A closed volume can be calculated as:
        // 1/3*surfIntegral(S*n dS)
        scalarField ownVolHole =
        -(1/3.0)*(((cf + Df) & nf * magSf));

        scalarField nbrVolHole =
        -(1/3.0)*(((nbrCf + nbrDf) & nbrNf * nbrMagSf));

        //- The volume is calculated as the difference between the volume inside
        //  the current patch minus the volume inside the nbr patch.
        scalarField volHole =  max(nbrVolHole + ownVolHole, 0.0);

        //- As the previous volume is calculated on both sides of the gap
        //  we need to divide by 2 when summing
        volume += 0.5*gSum(volHole);
        area += gSum(magSf);

        temperature += gSum(T*magSf);
        VoverT += 0.5*gSum(1/T*volHole);
    }

    //- Remove plena volumes
    gapV_ = volume;
    gapT_ = temperature/max(area, VSMALL);
    gapVoverT_ = VoverT;
}

void Foam::gapTRISO::correctBuffer()
{
    //- Initialize useful quantities
    scalar temperature(0.0);

    word bufferOuterName(gapGasDict_.lookupOrDefault<word>("bufferOuterPatches", "bufferOuter"));
    word bufferName(gapGasDict_.lookupOrDefault<word>("bufferZones", "Buffer"));
    label bufferOuterID(mesh_.boundaryMesh().findPatchID(bufferOuterName));
    label bufferID(mesh_.cellZones().findZoneID(bufferName));

    const scalarField T = T_.internalField();

    //- Beacuse the gas is inside the Buffer, so the volume-average temperature is used.
    const labelList& addrBuffer(mesh_.cellZones()[bufferID]);
    scalarField volumeBuffer(mesh_.V(), addrBuffer);
    scalarField temperatureBuffer(T, addrBuffer);
    temperature = gSum(temperatureBuffer*volumeBuffer);

    const fvPatch& p = mesh_.boundary()[bufferOuterID];
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(p);

    //- Patch center face vectors, displacement and temperature fields
    const vectorField& cf = patch.Cf();
    const fvPatchVectorField& Df = DorDD_.boundaryField()[bufferOuterID];
    const vectorField& Sf = patch.Sf();

    //- A closed volume can be calculated as:
    // 1/3*surfIntegral(S*n dS)
    scalarField VBufferKernel = mag(1/3.0*((cf + Df) & Sf));

    // for the inner faces collected in facezone
    const surfaceVectorField& meshCf(mesh_.Cf());
    const surfaceVectorField& meshSf(mesh_.Sf());
    const surfaceVectorField meshDf(fvc::interpolate(DorDD_));
    scalar VKernel = 0.0;



    forAll(mesh_.faceZones()[0], i)
    {

      label faceID(mesh_.faceZones()[0][i]);
      VKernel += mag(1/3.0*((meshCf[faceID] + meshDf[faceID]) & meshSf[faceID]));

    }

    //- Handle parallelisation
    reduce(VKernel, sumOp<scalar>());

    bufferV_ = max(gSum(VBufferKernel) - VKernel, 0.0) * porosity_;
    bufferT_ = temperature/max(gSum(volumeBuffer), SMALL);
    bufferVoverT_ = bufferV_ / bufferT_;
}


//- Override the same functions in gapGasModel. Beacuse in the case of TRISO, the
//- amount of gas could be very small at the beginning. So here instead of using
//- "SMALL" to verify the quantities, the "0.0" is used here.
void Foam::gapTRISO::correctMassFractions()
{
    scalar Ytot = sum(Y_);

    if (Ytot < 0.0) //modified part
    {
        FatalErrorIn
        (
            "Foam::gapFrapcon::correctMassFractions"
        ) << "Sum of mass fractions is zero or negative"
          << exit(FatalError);
    }
    else if (Ytot == 0.0)
    {
      forAll(Y_, i)
      {
          Y_[i] = 0.0;
      }
    }
    else
    {
      forAll(Y_, i)
      {
          Y_[i] /= Ytot;
      }
    }
}

void Foam::gapTRISO::correctMolFractions()
{
    scalar molTot(0.0);

    forAll(M_, i)
    {
        M_[i] = Y_[i]*gasM_/speciesW_[i];
        molTot += M_[i];
    }

    if (molTot < 0.0) //modified part
    {
        FatalErrorIn
        (
            "Foam::gapFrapcon::correctMolFractions"
        ) << "Sum of molar fractions is zero or negative"
          << exit(FatalError);
    }
    else if (molTot == 0.0)
    {
      M_ = 0.0;
    }
    else
    {
      M_ = M_/molTot;
    }
}


void Foam::gapTRISO::correctMass()
{
    //- Obtain list of gas components included in the fgr model
    const wordList& gasNames(fgr_.gasComponents());

    //- Obtain mol of gas released as calculated by fgr model
    const scalarList& deltaGasMol(fgr_.gasMols());

    //- Rescale mass fraction list using old time values
    Y_ = Y0_*gasM0_;

    //- Rescale current gas mass using old time value
    gasM_ = gasM0_;

    //- Loop over all gas species
    forAll(gasNames, gasNameI)
    {
            const word gasName = gasNames[gasNameI];
            bool found(false);

            forAll(species_, specieI)
            {
                    if(species_[specieI] == gasName)
                    {
                            Y_[specieI] +=
                            deltaGasMol[gasNameI]*speciesW_[specieI]/1000;

                            gasM_ +=
                            deltaGasMol[gasNameI]*speciesW_[specieI]/1000;

                            found = true;
                            break;
                    }
            }

            if(!found)
            {
                FatalErrorIn("gapGasModel::correctMass")
                << "Gas component "
                << gasNames[gasNameI] << nl
                << " not found in gapFrapcon" << abort(FatalError);
            }
    }

    // CO
    // After the update of fission gas, the update of CO is conducted here.
    label index = findIndex(species_, "CO");
    gasM_ -= Y_[index];
    Y_[index] = Mco_*speciesW_[index]/1000;
    gasM_ += Y_[index];

    correctMassFractions();
    correctMolFractions();
}


void Foam::gapTRISO::productionCO()
{
  scalar COprod; //The production of CO (atoms/fission)

  const volScalarField& heatSource = mesh_.lookupObject<volScalarField>(heatSourceName_);
  scalar total_fission_rate; // fisssion/s
  scalar total_CO; //The total CO (moles) till current time

  scalarField Q(heatSource);

  //- Do we need do the correction of volume ?
  total_fission_rate = gSum(Q*mesh_.V()/3.2e-11);

  if (mesh_.time().value() > currentTime1_ && total_fission_rate != 0.0)
  {
    currentTime1_ = mesh_.time().value();
    total_fission_ = total_fission_ + total_fission_rate * mesh_.time().deltaT().value();
  }

  if (CO_method_ == "Proksch")
  {
    COprod = ProkschModel();
  }
  else if (CO_method_ == "GA")
  {
    COprod = GAmodel();
  }
  else
  {
    FatalErrorIn("gapGasModel::productionCO()")
    << "The available CO production model are " << nl
    << "{" << nl << " Proksch" << nl << " GA" << nl << "}" << nl
    << "Please choose the correct one" << exit(FatalError);
  }




  total_CO = total_fission_ * COprod / 6.02214076e23;

  Mco_ = total_CO;
}

Foam::scalar Foam::gapTRISO::GAmodel()
{
  const surfaceScalarField& meshSf(mesh_.magSf());
  const surfaceScalarField& TKernelOuter(fvc::interpolate(T_));

  scalar temperature = 0.0;
  scalar area = 0.0;
  scalar Tav; // Surface average temperature of Kernel outer patch
  scalar COprod; //The production of CO (atoms/fission)
  //- Calculated the average surfacc temperatur of Kernel
  forAll(mesh_.faceZones()[0], i)
  {
    label faceID(mesh_.faceZones()[0][i]);
    temperature += TKernelOuter[faceID] * meshSf[faceID];
    area += meshSf[faceID];
  }
  //- Handle parallelisation
  reduce(temperature, sumOp<scalar>());
  reduce(area, sumOp<scalar>());
  Tav = temperature / area;
  COprod = min(0.625,  1.64*exp(-3311.0/Tav));

  return COprod;
}

Foam::scalar Foam::gapTRISO::ProkschModel()
{
  //Proksch model needs the TRISO particle surface temperature
  word OuterName(gapGasDict_.lookupOrDefault<word>("OuterPatch", "outer"));
  label OuterID(mesh_.boundaryMesh().findPatchID(OuterName));
  const fvPatchScalarField& T = T_.boundaryField()[OuterID];
  const fvPatch& p = mesh_.boundary()[OuterID];
  scalar Tav; // Surface average temperature of Kernel outer patch
  scalar COprod; //The production of CO (atoms/fission)
  scalar TintAv; // The time averaged surface temperature

  const scalarField& Sf = p.magSf();

  scalar temperature = gSum(T*Sf);
  scalar area = gSum(Sf);

  Tav = temperature / area;

  //static scalar Tint = Tav; //time integrated surface temperature

  //- Only do it when the code goes to the next step
  if (mesh_.time().value() > currentTime2_)
  {
    currentTime2_ = mesh_.time().value();
    Tint_ = Tint_ + Tav * mesh_.time().deltaT().value();
  }

  if (Tint_ < SMALL)
  {
    return 0.0;
  }

  TintAv = Tint_ / max(currentTime2_, SMALL);

  COprod = min(0.4,  pow(currentTime2_, 2)/(1.211e10*pow(10.0, 8500.0/TintAv)));

  return COprod;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gapTRISO::gapTRISO
(
    const fvMesh& mesh,
    const materials& mat,
    const volScalarField& T,
    const volVectorField& DorDD,
    const fissionGasRelease& fgr,
    const dictionary& gapGasDict
)
:
    gapGasModel(mesh, mat, T, DorDD, fgr, gapGasDict),
    heatSourceName_(gapGasDict.lookupOrDefault<word>("heatSourceName", "Q" )),
    CO_method_(gapGasDict.lookupOrDefault<word>("COProductionModel", "Proksch")),
    is3D_(gapGasDict.lookupOrDefault<bool>("model3D", false)),
    crackV_(0.0),
    gapOffsetV_(readScalar(gapGasDict.lookup("gapVolumeOffset"))),
    reserveV_(readScalar(gapGasDict.lookup("gasReserveVolume"))),
    porosity_(readScalar(gapGasDict.lookup("bufferPorosity"))),
    crackT_(290),
    reserveT_(290),
    bufferT_(290),
    reserveTlist_(),
    gapVoverT_(0),
    crackVoverT_(0),
    bufferVoverT_(0),
    currentTime1_(mesh_.time().value()),
    currentTime2_(mesh_.time().value()),
    total_fission_(0),
    total_fission0_(0),
    Tint_(0),
    Tint0_(0),
    Mco_(0.0)
{
  //***************************************************************************
  // Creat a faceZone between Buffer and Kernel;
  fvMesh& meshRef =  const_cast<fvMesh&>(mesh_);

  //- Find zone interfaces
  labelList zoneIDs(meshRef.nCells(), -1);
  //- Initialize index for faceZones
  label nFaceZones(meshRef.faceZones().size());

  forAll(meshRef.cellZones(), zoneI)
  {
      IndirectList<label>(zoneIDs, meshRef.cellZones()[zoneI]) = zoneI;
  }

  const labelList& owner = meshRef.owner();
  const labelList& neighbour = meshRef.neighbour();

  word bufferName(gapGasDict_.lookupOrDefault<word>("bufferZones", "Buffer"));
  label bufferID(mesh_.cellZones().findZoneID(bufferName));
  word kernelName(gapGasDict_.lookupOrDefault<word>("kernelZones", "Kernel"));
  label kernelID(mesh_.cellZones().findZoneID(kernelName));


  labelList interface;

  forAll(owner, faceI)
  {
    if (zoneIDs[owner[faceI]] == kernelID || zoneIDs[owner[faceI]] == bufferID)
    {
      if (zoneIDs[owner[faceI]] != zoneIDs[neighbour[faceI]])
      {
            interface.append(faceI);
      }
    }
  }
  boolList flipList(interface.size(), true);

  meshRef.faceZones().setSize(nFaceZones+1);
  meshRef.faceZones().set
  (
    nFaceZones,
    new faceZone
    (
      "bufferInner",
      interface,
      flipList,
      nFaceZones,
      meshRef.faceZones()
    )
  );
  //***************************************************************************
    if(gapGasDict.found("gasReserveTemperatureList"))
    {
#ifdef OPENFOAMFOUNDATION        
        reserveTlist_.set
        ( 
            new Function1s::Table<scalar>
            (
                "gasReserveTemperatureList", 
                gapGasDict.subDict("gasReserveTemperatureList")
            )
        );
#elif OPENFOAMESI        
        reserveTlist_.reset
        ( 
            Function1<scalar>::New
            (
                "gasReserveTemperatureList", 
                gapGasDict.subDict("gasReserveTemperatureList")
            )
        );
#endif        

        scalar currentTime(mesh_.time().value());

        reserveT_ =
        reserveTlist_->value(mesh_.time().timeToUserTime(currentTime));
    }
    else
    {
        reserveT_ = readScalar(gapGasDict.lookup("gasReserveTemperature"));
    }

    findPatchIDs("gap");

    //- Set the correction factor list. One list per patch (initial value is 1)
    scalingFactors_.setSize(mesh_.boundaryMesh().size());
    forAll(scalingFactors_, patchI)
    {
        scalingFactors_[patchI].setSize(mesh_.boundaryMesh()[patchI].size());
        scalingFactors_[patchI] = 1;
    }

    //- Derive species names from gapGas dictionary
    dictionary gapGas(readStream(this->typeName));
    close();

    const dictionary& dictMassFr = gapGas.subDict("massFractions");
    forAll(dictMassFr.toc(), entryI)
    {
        species_.append(dictMassFr.toc()[entryI]);
    }

    //- Read (or build default) dictionary for molar masses
    dictionary dictW("dictW");
    if(gapGasDict.found("speciesW"))
    {
        dictW = gapGasDict.subDict("speciesW");
    }
    else
    {
        dictW.add("Xe", 131.3);
        dictW.add("Ne", 20.183);
        dictW.add("Ar", 39.948);
        dictW.add("Kr", 83.8);
        dictW.add("Rn", 222);
        dictW.add("He", 4.0026);
        dictW.add("CO", 28.01);
    }

    //- Set species molar weight
    speciesW_.setSize(species_.size());

    forAll(species_, i)
    {
        if(dictW.found(species_[i]))
        {
            speciesW_[i] = readScalar(dictW.lookup(species_[i]));
        }
        else
        {
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl
            << " not found in speciesW in gapGasProperties" << exit(FatalError);
        }
    }

    //- Read (or build default) dictionary for conductiviy coefficient A
    dictionary dictAs("coefficientA");
    if(gapGasDict.found("conductivity_A"))
    {
        dictAs = gapGasDict.subDict("conductivity_A");
    }
    else
    {
        dictAs.add("Xe", 9.825e-5);
        dictAs.add("Ne", 5.000e-5);
        dictAs.add("Ar", 4.092e-4);
        dictAs.add("Kr", 1.966e-4);
        dictAs.add("Rn", 1.0);
        dictAs.add("He", 2.531e-3);
        dictAs.add("CO", 1.403e-4);
    }

    //- Set list of A coefficients for conductivity formula
    conductivityAs_.setSize(species_.size());
    forAll(species_, i)
    {
        if(dictAs.found(species_[i]))
        {
            conductivityAs_[i] = readScalar(dictAs.lookup(species_[i]));
        }
        else
        {
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl
            << " not found in conductivity_A in gapGasProperties" << exit(FatalError);
        }
    }

    //- Read (or build default) dictionary for conductiviy coefficient B
    dictionary dictBs("conductivity_B");
    if(gapGasDict.found("conductivity_B"))
    {
        dictBs = gapGasDict.subDict("conductivity_B");
    }
    else
    {
        dictBs.add("Xe", 0.7334);
        dictBs.add("Ne", 0.9228);
        dictBs.add("Ar", 0.6748);
        dictBs.add("Kr", 0.7006);
        dictBs.add("Rn", 1.0);
        dictBs.add("He", 0.7146);
        dictBs.add("CO", 0.9090);
    }

    //- Set list of B coefficients for conductivity formula
    conductivityBs_.setSize(species_.size());
    forAll(species_, i)
    {
        if(dictAs.found(species_[i]))
        {
            conductivityBs_[i] = readScalar(dictBs.lookup(species_[i]));
        }
        else
        {
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl
            << " not found in conductivity_B in gapGasProperties" << exit(FatalError);
        }
    }

    //- Read current mass, mass fractions, pressure
    //  and calculate molar fractions and gap gas volumes
    readData(readStream(this->typeName));
    close();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gapTRISO::~gapTRISO()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::gapTRISO::kappa(const scalar T) const
{
    scalarList ks_(species_.size(), 0.0);
    scalar Ytot = sum(Y_);
    //In the case no gas is presence, the thermal conductivity is 0.0;
    if (Ytot == 0.0)
    {
      return 0.0;
    }

    ks_ = conductivityAs_*pow(T, conductivityBs_);

    scalar kMix(0.0);

    forAll(speciesW_, specieI)
    {
            scalarList deltas_(species_.size(), 0.0);
            deltas_[specieI] = 1.0;
            scalar summTerm(0.0);

            scalar kI = ks_[specieI];
            scalar nI = M_[specieI];

            forAll(species_, specieJ)
            {
                    scalar phiIJ(0.0);
                    scalar psiIJ(0.0);
                    scalar kJ = ks_[specieJ];
                    scalar MI = speciesW_[specieI]/1000;
                    scalar MJ = speciesW_[specieJ]/1000;

                    phiIJ =
                    (
                        pow( (1 + pow(kI/kJ, 0.5)*pow(MI/MJ, 0.25)), 2.0)
                    )/(pow(2, 1.5)*pow(1 + MI/MJ , 0.5));

                    psiIJ = phiIJ*
                    (
                        1 + 2.41*
                        (
                            (MI - MJ)*(MI - 0.142*MJ)
                        )/(pow(MI + MJ, 2.0))
                    );

                    summTerm += (1-deltas_[specieJ])*psiIJ*M_[specieJ];

            }

            kMix += ( (kI*nI)/(nI + summTerm) ) ;
    }

    return kMix;
}


Foam::scalar Foam::gapTRISO::a(const scalar T) const
{
    scalar a_He = 0.425 - (2.3e-4*min(T, 1300));
    scalar a_Xe = 0.749 - (2.5e-4*min(T, 1300));
    const scalar W_He = 4.0026;
    const scalar W_Xe = 131.3;

    scalar a_mix = 0;

    forAll(species_, i)
    {
        scalar W_i = speciesW_[i];
        scalar a_i = a_He + (W_i - W_He)/(W_Xe - W_He)*(a_Xe - a_He);

        a_mix += M_[i]*a_i/sqrt(W_i);
    }

    return a_mix;
}


void Foam::gapTRISO::correct()
{
    //- Check fgr model for additional gas release and update gas mass,
    //  mass fractions and molar fractions
    productionCO();
    correctMass();
    correctGap();
    correctBuffer();


    //- Update reserve temperature (if time list was given)
    if(reserveTlist_.valid())
    {
        scalar currentTime(mesh_.time().value());

        reserveT_ =
        reserveTlist_->value(mesh_.time().timeToUserTime(currentTime));
    }

    //- Update pressure
    if(gasPType_ == "fromModel")
    {
        //- Total number of moles
        scalar molN(0);

        forAll(Y_, nuclideI)
        {
            molN += Y_[nuclideI]*gasM_/(speciesW_[nuclideI]/1000);
        }

        //- Universal gas constant
        scalar R = Foam::constant::physicoChemical::R.value();

        //- Update pressure (ideal gas)
        gasP_ = molN*R/(gapVoverT_ + gapOffsetV_/max(gapT_,VSMALL)
                      + bufferVoverT_ + reserveV_/reserveT_ + crackVoverT_);
    }
    else if(gasPType_ == "fromList")
    {
        scalar currentTime(mesh_.time().value());

        gasP_ =
        gasPressureList_->value(mesh_.time().timeToUserTime(currentTime));
    }
    else if(gasPType_ == "fixed")
    {
        //gapP does not change
    }
    else
    {
        FatalErrorIn("gapTRISO:correct()")
        << "gasPressureType: "
        << gasPType_
        << " not known. Use either:" << nl
        << "- \"fromModel\"" << nl
        << "- \"fromList\"" << nl
        << "- \"fixed\"" << exit(FatalError);
    }

    //*******Output relevant information *****/

    Info<< "Gap conditions:" << endl
        << tab << "Gap volume: " << gapV_ << endl
        << tab << "Gap temperature: " << gapT_ << endl
        << tab << "Buffer gas volume: " << bufferV_ << endl
        << tab << "Buffer gas average temperature: " << bufferT_ << endl
        << tab << "Reserve volume: " << reserveV_ << endl
        << tab << "Reserve temperature: " << reserveT_ << endl
        << tab << "Crack gas volume: " << crackV_ << endl
        << tab << "Crack gas average temperature: " << crackT_ << endl
        << tab << "Gas pressure: " << gasP_ << endl
        << tab << "Gas mass: " << gasM_ << endl
        << tab << "Y: " << Y_ << endl
        << tab << "M: " << M_ << endl

        << endl;
}


bool Foam::gapTRISO::readData(Istream& is)
{
    dictionary dict(is);

    // The gas mass is either calculated or given in gapGas dict.
    // V, T and p initial conditions are necessary for calculating
    // initial gass mass and number of mols.

    // Read pressure
    gasPType_ = dict.lookupOrDefault<word>("gasPressureType", "fromModel");
    if(gasPType_ == "fromModel" or gasPType_ == "fixed")
    {
        gasP_ = readScalar(dict.lookup("gasPressure"));
    }
    else if(gasPType_ == "fromList")
    {
#ifdef OPENFOAMFOUNDATION        
        gasPressureList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "gapPressureList", 
                dict.subDict("gapPressureList")
            )
        );
#elif OPENFOAMESI        
        gasPressureList_.reset
        ( 
            Function1<scalar>::New
            (
                "gapPressureList", 
                dict.subDict("gapPressureList")
            )
        );
#endif           

        scalar currentTime(mesh_.time().value());

        gasP_ =
        gasPressureList_->value(mesh_.time().timeToUserTime(currentTime));
    }
    else
    {
        FatalErrorIn("gapTRISO:readData")
        << "Gas pressure type "
        << gasPType_ << nl
        << " not known. Use fromModel, fromList or fixed" << exit(FatalError);
    }


    // Gap volume, temperature and V/T contribution
    correctGap();
    correctBuffer();
    // FatalErrorIn("gapTRISO:readData")
    // << "Here" << exit(FatalError);

//- Read mass fraction list
    const dictionary& dictMassFr = dict.subDict("massFractions");

    Y_.resize(species_.size());

    forAll(species_, i)
    {
        if(dictMassFr.found(species_[i]))
        {
            Y_[i] = readScalar(dictMassFr.lookup(species_[i]));
        }
        else
        {
            FatalErrorIn("gapGasModel")
            << "Gas component "
            << species_[i] << nl
            << " not found in massFractions in gapGas" << exit(FatalError);
        }
    }

    correctMassFractions();

    Y0_ = Y_;

    //- Read mass
    gasM0_ = dict.lookupOrDefault<scalar>("gasMass", 0.0);
    //- if initial mass is not present, calculate it
    if( gasM0_ <= 0.0 )
    {
            calcInitialMass();
    }

    gasM_ = gasM0_;

    //- Set molar fraction list
    M_.resize(species_.size());
    correctMolFractions();
    M0_ = M_;

    //- Read total fission
    total_fission0_ = dict.lookupOrDefault<scalar>("totalFission", 0.0);
    total_fission_ = total_fission0_;

    Tint0_ = dict.lookupOrDefault<scalar>("integratedTemperature", 0.0);
    Tint_ = Tint0_;

    Info<< "Gap conditions:" << endl
        << tab << "Gap volume: " << gapV_ << endl
        << tab << "Gap temperature: " << gapT_ << endl
        << tab << "Buffer gas volume: " << bufferV_ << endl
        << tab << "Buffer gas average temperature: " << bufferT_ << endl
        << tab << "Reserve volume: " << reserveV_ << endl
        << tab << "Reserve temperature: " << reserveT_ << endl
        << tab << "Crack gas volume: " << crackV_ << endl
        << tab << "Crack gas average temperature: " << crackT_ << endl
        << tab << "Gas pressure: " << gasP_ << endl
        << tab << "Gas mass: " << gasM_ << endl
        << tab << "Y: " << Y_ << endl
        << tab << "M: " << M_ << endl
        << endl;

    return !is.bad();
}


bool Foam::gapTRISO::writeData(Ostream& os) const
{
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    // Obtain model angular fraction to calculate the gap info for the
    // corresponding  3-D sphere
    scalar fraction = 1.0;
    if (is3D_ == false)
    {
      scalar angularFraction = globalOpt.angularFraction();
      fraction = pow(degToRad(angularFraction*360),2)/4/M_PI;
    }
    else
    {
      // TO DO: other cases, e.x, a quarter.
      fraction = 1.0/8.0;
    }

    dictionary data;

    forAll(species_, i)
    {
        data.add(species_[i], Y_[i]);
    }

    os.writeKeyword("massFractions") << data << token::END_STATEMENT
    << nl << nl;

    os.writeKeyword("gasPressureType") << gasPType_ << token::END_STATEMENT
    << nl << nl;

    os.writeKeyword("gapVolume (in the OFFBEAT model)")
    << gapV_ << token::END_STATEMENT << nl;
    os.writeKeyword("gapVolume (in the corresponding whole sphere)")
    << gapV_ / fraction << token::END_STATEMENT << nl;
    os.writeKeyword("gapTemperature") << gapT_ << token::END_STATEMENT
    << nl << nl;

    os.writeKeyword("BufferGasVolume (in the OFFBEAT model)")
    << bufferV_ << token::END_STATEMENT << nl;
    os.writeKeyword("BufferGasVolume (in the corresponding whole sphere)")
    << bufferV_ / fraction << token::END_STATEMENT << nl;
    os.writeKeyword("BufferGasTemperature") << bufferT_ << token::END_STATEMENT
    << nl << nl;

    os.writeKeyword("reserveVolume (in the OFFBEAT model)")
    << reserveV_ << token::END_STATEMENT << nl;
    os.writeKeyword("reserveVolume (in the corresponding whole sphere)")
    << reserveV_/fraction << token::END_STATEMENT << nl;
    os.writeKeyword("reserveT_") << reserveT_ << token::END_STATEMENT
    << nl << nl;

    os.writeKeyword("crackGasVolume (in the OFFBEAT model)")
    << crackV_ << token::END_STATEMENT << nl;
    os.writeKeyword("crackGasVolume (in the corresponding whole sphere)")
    << crackV_/fraction << token::END_STATEMENT << nl;
    os.writeKeyword("crackGasTemperature") << crackT_ << token::END_STATEMENT
    << nl << nl;

    os.writeKeyword("gasPressure") << gasP_ << token::END_STATEMENT << nl;
    os.writeKeyword("gasMass") << gasM_ << token::END_STATEMENT << nl;
    os.writeKeyword("totalFission") << total_fission_ << token::END_STATEMENT << nl;
    os.writeKeyword("integratedTemperature") << Tint_ << token::END_STATEMENT << nl;

    return os.good();
}

// ************************************************************************* //
