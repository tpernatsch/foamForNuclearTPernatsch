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

#include "fissionProductsDiffusionSolver.H"
#include "zeroCurrentActinidesRedistributionFvPatchScalarField.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fissionProductsDiffusionSolver, 0);
    addToRunTimeSelectionTable
    (
        fissionProductsDiffusionSolver,
        fissionProductsDiffusionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fissionProductsDiffusionSolver::initializeFpConcentration(const dictionary& FpOptDict)
{
    const labelListList& matAddrList(mat_.matAddrList());
    const dictionary& Fp = FpOptDict.subDict("intialFpConcentration");

    if (mesh_.time().value() == 0) //for restart
    {
      //TO DO? : For now only Kernel can have the initial concentration. Maybe
      //we should do the same thing for other parts?
      forAll(mat_.materialsList(), i)
      {
        if (mat_.materialsList()[i].name() == "UO2")
        {
          const labelList& addr(matAddrList[i]);
          const fuelMaterial& fuelMat = refCast<const fuelMaterial>(mat_.materialsList()[i]);
          enrichment_ = fuelMat.enrichment();

          forAll(addr, j)
          {
              const label cellI = addr[j];

              forAll(FpNames_, k)
              {
                Fp_[k]()[cellI] = readScalar(Fp.lookup(FpNames_[k]));
              }

          }
       }
      }
    }
    //Calculate Initial Fission products quantities
    scalarField volumeMat(mesh_.V());

    forAll(FpNames_, k)
    {
      totalInital_[k] = gSum(volumeMat*Fp_[k]());
      Fp_[k]().correctBoundaryConditions();
    }
}

void Foam::fissionProductsDiffusionSolver::updateSourceTerm
(
  const labelList& addr,
  const scalarField& Bu,
  const scalarField& f,
  const scalar& e
)
{
  const scalar Na = 6.0221408e23;
  forAll(FpNames_, i)
  {
    if (FpNames_[i] == "Cs")
    {
      forAll(addr, j)
      {
        // atom/fission * fission/s/m³ / mol^-1 = mol/m³/s
        s_[i]()[addr[j]] = (e < 0.175 ? 0.14 : 0.16) * f[addr[j]] / Na;
      }
    }
    else if (FpNames_[i] == "Kr")
    {
      forAll(addr, j)
      {
        s_[i]()[addr[j]] = 0.297 * f[addr[j]] / Na; //Kr+Xe
      }
    }
    else if (FpNames_[i] == "Ag")
    {
      forAll(addr, j)
      {
        const scalar b = max(1.0, Bu[addr[j]]*1e-3/0.881/9.3706); //Convert Bu in MWd/tUO2 to %FIMA
        if (e < 0.175)
        {
          s_[i]()[addr[j]] = 1.31625e-3*pow(b, 0.55734) * f[addr[j]] / Na;
        }
        else
        {
          s_[i]()[addr[j]] = 8.24492e-4*pow(b, 0.53853) * f[addr[j]] / Na;
        }
      }
    }
  }
}

void Foam::fissionProductsDiffusionSolver::updateCoefficients()
{
    const scalarField& Bui(Bu_.internalField());
    const scalarField& Qi(Q_.internalField());

    // W/m³ / (J/fission) = fission/m³/s;
    const scalarField f = Qi / 3.2e-11;

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
      const labelList& addr(matAddrList[i]);

      if(mat_.materialsList()[i].name() == "UO2")
      {
         updateSourceTerm(addr, Bui, f, enrichment_);
      }

       forAll(FpNames_, j)
       {
         diffModels_[i] -> updateCoef(diff_[j](), T_, addr, FpNames_[j]);
       }
    }

    forAll(FpNames_, k)
    {
      diff_[k]().correctBoundaryConditions();
    }
}
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fissionProductsDiffusionSolver::fissionProductsDiffusionSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    transportSolver(mesh, mat, elementTransportDict, solverName),
    // regIOobject
    // (
    //     IOobject
    //     (
    //         "CsTransport",
    //         mesh.time().timeName()/"uniform",
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     )
    // ),
    T_(mesh_.lookupObject<volScalarField>("T")),
    phi_(mesh_.lookupObject<volScalarField>("fastFluence")),
    Bu_(mesh_.lookupObject<volScalarField>("Bu")),
    Q_(mesh_.lookupObject<volScalarField>("Q")),
    releasedInc_(),
    releasedPrev_(),
    released_(),
    currentTime_(),
    totalInital_(),
    totalProd_(),
    totalProdPre_(),
    patchName_(),
    FpNames_(),
    FpName_("FP")
{
  dictionary FpDict
  (
      elementTransportDict.
      subOrEmptyDict("FpDiffusionOptions")
  );

  patchName_ = FpDict.lookupOrDefault<word>("outerPatches", "outer");
  const dictionary& Fp = FpDict.subDict("intialFpConcentration");
//  List<autoPtr<volScalarField>>FpList(2);
  //FpList.resize(FpNames_.size());

  forAll(Fp.toc(), entryI)
  {
    FpNames_.append(Fp.toc()[entryI]);
  }

  //Decide the type of Fp to solve according to the input. So we can save the
  //the computational resources from the non-needs Fp.
  Fp_.resize(FpNames_.size());
  diff_.resize(FpNames_.size());
  s_.resize(FpNames_.size());

  forAll(FpNames_, i)
  {
    Fp_[i].set
    (
      new volScalarField
      (
        IOobject
        (
            FpNames_[i],
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimMoles/dimVolume, 0.0), // mol/m³
        "fixedValue"
      )
    );

    diff_[i].set
    (
      new volScalarField
      (
        IOobject
        (
            "diff"+FpNames_[i],
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE //change back
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, 0.0), // m²/s
        "calculated"
      )
    );

    s_[i].set
    (
      new volScalarField
      (
        IOobject
        (
            "source"+FpNames_[i],
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE //change back
        ),
        mesh_,
        dimensionedScalar("", dimMoles/dimVolume/dimTime, 0.0), // mol/m³/s
        "calculated"
      )
    );
  }

  releasedInc_.resize(FpNames_.size());
  releasedPrev_.resize(FpNames_.size());
  released_.resize(FpNames_.size());
  totalInital_.resize(FpNames_.size());
  totalProd_.resize(FpNames_.size());
  totalProdPre_.resize(FpNames_.size());
  Fpleft_.resize(FpNames_.size());


  diffModels_.resize(mat_.materialsList().size());


  forAll(mat_.materialsList(), i)
  {
    const materialModel& Mat = refCast<const materialModel>(mat_.materialsList()[i]);
    dictionary MatDict = Mat.materialModelDict();
    if (mat_.materialsList()[i].name() == "UO2")
    {
      diffModels_[i] = diffCoefModel::New(mesh_, MatDict, "kernelArrhenius");
    }
    else if (mat_.materialsList()[i].name() == "PyC")
    {
      diffModels_[i] = diffCoefModel::New(mesh_, MatDict, "PyCArrhenius");
    }
    else if (mat_.materialsList()[i].name() == "SiC")
    {
      diffModels_[i] = diffCoefModel::New(mesh_, MatDict, "SiCArrhenius");
    }
    else if (mat_.materialsList()[i].name() == "graphite")
    {
      diffModels_[i] = diffCoefModel::New(mesh_, MatDict, "graphiteArrhenius");
    }
    else
    {
      diffModels_[i] = diffCoefModel::New(mesh_, MatDict, "inputArrhenius");
    }
  }







    // Initialize Am concentration
    initializeFpConcentration(FpDict);


    // Correct coefficients
    updateCoefficients();




}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::fissionProductsDiffusionSolver::correctFp
(
  volScalarField& Fp,
  const volScalarField diffFp,
  const volScalarField s,
  const int nCorr
)
{
  int nInnerIter = 0;

  bool convergedInner(false);
  do
  {
      // Store prev iter values for rel residual
      Fp.storePrevIter();

      // Equation for Fp diffusion, the decay term is not yet considered.
      fvScalarMatrix FpEqn
      (
            fvm::ddt(Fp)
          - fvm::laplacian(diffFp, Fp, "laplacianFp")
          ==
          s
      );


      // Relax equation
      FpEqn.relax();

      // Compute residual
      residual_ = FpEqn.solve().max().initialResidual();

      // Store initial residual if in first inner iteration
      if(nInnerIter==0)
      {
          initialResidual_ = residual_;
      }

      // Relax fields
      Fp.relax();

      // Calculate a different residual based on the relative change
      scalar denom =
      gMax(mag(Fp.primitiveField() - Fp.oldTime().primitiveField()));

      scalar num =
      gMax(mag(Fp.primitiveField() - Fp.prevIter().internalField()));

      if (denom < SMALL)
      {
          denom = max(gMax(mag(Fp.primitiveField())), SMALL);
      }

      relResidual_ = num/denom;

      if(denom < 1e-5)
      {
          relResidual_ = VSMALL;
      }

      Info << "relResidual" << Fp.name() << " " << relResidual_ << endl;

      //We just use one criteria for all the fission products.
      convergedInner = converged(fieldName());

  } while
  (
      not(convergedInner)
      && ++nInnerIter < nCorr
  );
}

void Foam::fissionProductsDiffusionSolver::correct()
{
    // Update thermal diffusion coefficient
    updateCoefficients();
    scalarList residualFp(FpNames_.size());

    const dictionary& stressControl =  mesh_.solutionDict().subDict("stressAnalysis");
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    if (mesh_.time().value() > currentTime_)
    {
      releasedPrev_ = released_;
      totalProdPre_ = totalProd_;
      currentTime_ = mesh_.time().value();
    }


    forAll(FpNames_, i)
    {
      correctFp(Fp_[i](), diff_[i](), s_[i](), nCorr);
      residualFp[i] = initialResidual_;
    }

    initialResidual_ = max(residualFp);




    outputInfo();


    //Out put information

    Info<< "Fission products transport:" << endl;
    forAll(FpNames_, i)
    {
      Info << tab << "Initial " << FpNames_[i] << " amount: " << totalInital_[i] << endl
      << tab << "Production of " << FpNames_[i] << ": " << totalProd_[i] << endl
      << tab << "Released " << FpNames_[i] << " amount: " << released_[i] << endl
      << tab << "Remained " << FpNames_[i] << " amount: " << Fpleft_[i] << endl
      << tab << "Fraction of " << FpNames_[i] << " release: " << released_[i]/max(VSMALL,(totalProd_[i]+totalInital_[i])) << endl;
    }

    Info << endl;

}


void Foam::fissionProductsDiffusionSolver::outputInfo()
{
  label patchID(mesh_.boundaryMesh().findPatchID(patchName_));

  if(patchID == -1)
  {
      FatalErrorIn("Foam::fissionProductsDiffusionSolver::outputInfo()")
      << patchName_ << " Patch"<< " not found in mesh." << exit(FatalError);
  }
  const fvPatch& p = mesh_.boundary()[patchID];
  const vectorField& Sf = p.Sf();
  scalarField volumeMat(mesh_.V());
  forAll(FpNames_, i)
  {
    const volVectorField flux = fvc::grad(Fp_[i]())*diff_[i]();
    const fvPatchVectorField& fluxf = flux.boundaryField()[patchID];
    releasedInc_[i] = gSum(-fluxf & Sf);
    released_[i] = releasedPrev_[i] + releasedInc_[i]*mesh_.time().deltaT().value();

    scalarField volume(mesh_.V());

    totalProd_[i] = totalProdPre_[i] + gSum(volume*s_[i]())*mesh_.time().deltaT().value();
    Fpleft_[i] = gSum(volumeMat*Fp_[i]());
  }

}

// bool Foam::fissionProductsDiffusionSolver::writeData(Ostream& os) const
// {
//   Info << " writeData" << endl;
//   os.writeKeyword("Realised fission products") << released_ << token::END_STATEMENT<< nl << nl;
//   os.writeKeyword("Initial fission products quantities") << totalInital_ << token::END_STATEMENT<< nl << nl;
//
//   return os.good();
// }

// ************************************************************************* //
