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

#include "gapGasTimeTabulated.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gapGasTimeTabulated, 0);
    addToRunTimeSelectionTable
    (
        gapGasModel, 
        gapGasTimeTabulated,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::gapGasTimeTabulated::correctMassFractions()
{
    scalar Ytot = sum(Y_);
    if (Ytot < SMALL)
    {
        FatalErrorIn
        (
            "Foam::gapGasTimeTabulated::correctMassFractions"
        ) << "Sum of mass fractions is zero or negative"
          << exit(FatalError);
    }

    forAll(Y_, i)
    {
        Y_[i] /= Ytot;
    }
}


void Foam::gapGasTimeTabulated::correctMolFractions()
{
    scalar molTot(0.0);

    forAll(M_, i)
    {
        M_[i] = Y_[i]*gasM_/speciesW_[i];
        molTot += M_[i];
    }

    if (molTot < SMALL)
    {
        FatalErrorIn
        (
            "Foam::gapGasTimeTabulated::correctMolFractions"
        ) << "Sum of molar fractions is zero or negative"
          << exit(FatalError);
    }

    M_ = M_/molTot;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gapGasTimeTabulated::gapGasTimeTabulated
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
    pressureData_(gapGasDict.lookup("gapPressures")),
    ArgonData_(gapGasDict.lookup("Ar")),
    HeliumData_(gapGasDict.lookup("He")),
    KryptonData_(gapGasDict.lookup("Kr")),
    NeonData_(gapGasDict.lookup("Ne")),
    RadonData_(gapGasDict.lookup("Rn")),
    XenonData_(gapGasDict.lookup("Xe")),
    timeData_(gapGasDict.lookup("timePoints")),
    timeMethod_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            gapGasDict.lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    boundsMethod_
    (
        interpolateTableBase::outOfBoundsMethodNames_
        [
            gapGasDict.lookupOrDefault<word>
            ("outOfBoundsMethod", "error")
        ]
    ),
    pressureTable_(timeData_, pressureData_, timeMethod_, boundsMethod_),
    ArgonTable_(timeData_, ArgonData_, timeMethod_, boundsMethod_),
    HeliumTable_(timeData_, HeliumData_, timeMethod_, boundsMethod_),
    KryptonTable_(timeData_, KryptonData_, timeMethod_, boundsMethod_),
    NeonTable_(timeData_, NeonData_, timeMethod_, boundsMethod_),
    RadonTable_(timeData_, RadonData_, timeMethod_, boundsMethod_),
    XenonTable_(timeData_, XenonData_, timeMethod_, boundsMethod_)
{
    // Derive species names from gapGas dictionary
    dictionary gapGas(readStream(this->typeName));
    close();

    const dictionary& dictMassFr = gapGas.subDict("massFractions");
    forAll(dictMassFr.toc(), entryI)
    {
        species_.append(dictMassFr.toc()[entryI]);
    }

    // Read (or build default) dictionary for molar masses
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
    }

    // Set species molar weight
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

    // Read (or build default) dictionary for conductiviy coefficient A
    dictionary dictAs("coefficientA");
    if(gapGasDict.found("conductivity_A"))
    {
        dictAs = gapGasDict.subDict("conductivity_A");
    }
    else
    {
        dictAs.add("Xe", 9.825e-5);
        dictAs.add("Ne", 1.0);
        dictAs.add("Ar", 4.092e-4);
        dictAs.add("Kr", 1.966e-4);
        dictAs.add("Rn", 1.0);
        dictAs.add("He", 2.531e-3);
    }

    // Set list of A coefficients for conductivity formula
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

    // Read (or build default) dictionary for conductiviy coefficient B
    dictionary dictBs("conductivity_B");
    if(gapGasDict.found("conductivity_B"))
    {
        dictBs = gapGasDict.subDict("conductivity_B");
    }
    else
    {
        dictBs.add("Xe", 0.7334);
        dictBs.add("Ne", 1.0);
        dictBs.add("Ar", 0.6748);
        dictBs.add("Kr", 0.7006);
        dictBs.add("Rn", 1.0);
        dictBs.add("He", 0.7146);
    }

    // Set list of B coefficients for conductivity formula
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

    // Read current mass, mass fractions, pressure
    // and calculate molar fractions and gap gas volumes
    readData(readStream(this->typeName));
    close();

    // Read the gas presure and mass fraction histories
    readHistories();

    // TODO
    // Temporary fix to change the write path of the file gapGas in time/uniform
    // folders even though it is initialized in the 0/ folder.
    // In the future it will be inizialized in 0/uniform.
    fileName& localPath = const_cast<fileName&>(this->local());
    localPath = "uniform";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gapGasTimeTabulated::~gapGasTimeTabulated()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::scalar Foam::gapGasTimeTabulated::kappa(const scalar T) const
{
    scalarList ks_(species_.size(), 0.0);

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

Foam::scalar Foam::gapGasTimeTabulated::a(const scalar T) const
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

void Foam::gapGasTimeTabulated::correct()
{  
    correctMassFractions();
    correctMolFractions();

    // Update pressure
    scalar currentTime(mesh_.time().value());

    //gasP_ = gasPressureList->value(mesh_.time().timeToUserTime(currentTime));
    gasP_ = lookupPressure(currentTime);

    // Update mass fraction list
/*
    Ar_ = ArgonList->value(mesh_.time().timeToUserTime(currentTime));
    He_ = HeliumList->value(mesh_.time().timeToUserTime(currentTime));
    Kr_ = KryptonList->value(mesh_.time().timeToUserTime(currentTime));
    Ne_ = NeonList->value(mesh_.time().timeToUserTime(currentTime));
    Rn_ = RadonList->value(mesh_.time().timeToUserTime(currentTime));
    Xe_ = XenonList->value(mesh_.time().timeToUserTime(currentTime));
*/
    Ar_ = lookupAr(currentTime);
    He_ = lookupHe(currentTime);
    Kr_ = lookupKr(currentTime);
    Ne_ = lookupNe(currentTime);
    Rn_ = lookupRn(currentTime);
    Xe_ = lookupXe(currentTime);

    // Add mass fraction to species list
    Y_ = {Ar_, He_, Kr_, Ne_, Rn_, Xe_};

    // Output relevant information //

    Info<< "Gap conditions:" << endl
        << tab << "Gas pressure: " << gasP_ << endl
        << tab << "Y: " << Y_ << endl
        << tab << "M: " << M_ << endl
        << endl;
}

void Foam::gapGasTimeTabulated::readHistories()
{
    scalar currentTime(mesh_.time().value());
    // Read pressure
    gasP_ = lookupPressure(currentTime);

    // Read mass fraction list
    Ar_ = lookupAr(currentTime);
    He_ = lookupHe(currentTime);
    Kr_ = lookupKr(currentTime);
    Ne_ = lookupNe(currentTime);
    Rn_ = lookupRn(currentTime);
    Xe_ = lookupXe(currentTime);

    // Add mass fraction to species list
    Y_ = {Ar_, He_, Kr_, Ne_, Rn_, Xe_};

    correctMassFractions();
    Y0_ = Y_;

    // Read mass
    //gasM0_ = dict.lookupOrDefault<scalar>("gasMass", 1.0);
    gasM0_ = 1.0; //initial value; will be updated in the next iterations
    gasM_ = gasM0_;

    // Set molar fraction list
    M_.resize(species_.size());
    correctMolFractions();
    M0_ = M_;

    Info<< "Gap conditions:" << endl
        << tab << "Gas pressure: " << gasP_ << endl
        << tab << "Gas mass: " << gasM_ << endl
        << tab << "Y: " << Y_ << endl
        << tab << "M: " << M_ << endl
        << endl;

}

bool Foam::gapGasTimeTabulated::readData(Istream& is)
{
    return !is.bad();
}

bool Foam::gapGasTimeTabulated::writeData(Ostream& os) const
{    
    dictionary data;
    
    forAll(species_, i)
    {
        data.add(species_[i], Y_[i]);
    }
    
    os.writeKeyword("massFractions") << data << token::END_STATEMENT
    << nl << nl;
    os.writeKeyword("gasPressure") << gasP_ << token::END_STATEMENT << nl;

    return os.good();
}

// ************************************************************************* //
