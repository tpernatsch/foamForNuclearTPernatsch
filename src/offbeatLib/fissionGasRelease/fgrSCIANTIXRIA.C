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

#include "fgrSCIANTIXRIA.H"
#include "addToRunTimeSelectionTable.H"
#include "gapGasModel.H"
#include "zeroGradientFvPatchFields.H"
#include "globalFieldLists.H"
#include <string>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fgrSCIANTIXRIA, 0);
    addToRunTimeSelectionTable(fissionGasRelease, fgrSCIANTIXRIA, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fgrSCIANTIXRIA::fgrSCIANTIXRIA
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& fissionGasReleaseDict
)
:
    fissionGasRelease(mesh, mat, fissionGasReleaseDict),
    regIOobject
    (
        IOobject
        (
            "SCIANTIX_Fields",
            mesh.time().timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    sciantixFieldsDict_
    (
        headerOk() ? 
        readStream(this->typeName) : 
        dictionary::null
    ),
    addToRegistry_
    (       
#ifdef OPENFOAMFOUNDATION
        fissionGasReleaseDict_.isNull() ? 
#elif OPENFOAMESI
        fissionGasReleaseDict_.isNullDict() ? 
#endif             
        false :
        fissionGasReleaseDict_.lookupOrDefault("addToRegistry", false)
       
    ),
    Gas_produced_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_produced", sciantixFieldsDict_, addToRegistry_)
    ),
    Gas_grain_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_grain", sciantixFieldsDict_, addToRegistry_)),
    Gas_grain_solution_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_grain_solution", sciantixFieldsDict_, addToRegistry_)
    ),
    Gas_grain_bubbles_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_grain_bubbles", sciantixFieldsDict_, addToRegistry_)
    ),
    Gas_boundary_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_boundary", sciantixFieldsDict_, addToRegistry_)
    ),
    Gas_released_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_released", sciantixFieldsDict_, addToRegistry_)
    ),
    Gas_released_t0_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Gas_released_t0", sciantixFieldsDict_, addToRegistry_)
    ),
    Effective_burn_up_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Effective_burn_up", sciantixFieldsDict_, addToRegistry_)
    ),
    Helium_produced_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Helium_produced", sciantixFieldsDict_, addToRegistry_)
    ),
    Helium_grain_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Helium_grain", sciantixFieldsDict_, addToRegistry_)
    ),
    Helium_grain_solution_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Helium_grain_solution", sciantixFieldsDict_, addToRegistry_)
    ),
    Helium_grain_bubbles_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Helium_grain_bubbles", sciantixFieldsDict_, addToRegistry_)
    ),
    Helium_boundary_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Helium_boundary", sciantixFieldsDict_, addToRegistry_)
    ),
    Helium_released_
    (
        createOrRead<scalar>
        (mesh, "SCIANTIX::Helium_released", sciantixFieldsDict_, addToRegistry_)
    ),
    relax_
    (       
#ifdef OPENFOAMFOUNDATION
        fissionGasReleaseDict_.isNull() ? 
#elif OPENFOAMESI
        fissionGasReleaseDict_.isNullDict() ? 
#endif             
        1 :
        fissionGasReleaseDict_.lookupOrDefault("relax", 1.0)
    ),
    releaseHBS_(true),
    buReleaseThresholdHBS_(80000),
    temperatureReleaseThresholdHBS_(1000),
    damageReleaseThreshold_(0.85)
{ 
    // Calculate current fgr value
    scalarField produced((mesh_.V()*Gas_produced_*pos(grainRadius_)));
    scalarField released((mesh_.V()*Gas_released_*pos(grainRadius_)));

    // Gas released in RIA is calculated as for refabricated rods (i.e 
    // relative to initial inventory). Here we are setting 
    // Gas_released_t0_ = Gas_released_ for cases where Gas_released_t0_
    // has not been calculated
    if(gSum(Gas_released_t0_) < SMALL)
    {
        Gas_released_t0_ = Gas_released_;
    }

    // TODO: Using pos(grainRadius_) is a trick that works for now but it
    // is not safe
    scalarField released0((mesh_.V()*Gas_released_t0_*pos(grainRadius_)));
    fgr_ = (released - released0)/max(produced - released0, SMALL)*100;
    fgr0_ = fgr_;

    // Set release thresholds
    if 
    (
#ifdef OPENFOAMFOUNDATION
        fissionGasReleaseDict_.isNull()
#elif OPENFOAMESI
        fissionGasReleaseDict_.isNullDict()
#endif
    )
    {    
        releaseHBS_ = 
        fissionGasReleaseDict.lookupOrDefault("releaseHBS", true);
        buReleaseThresholdHBS_ = 
        fissionGasReleaseDict_.lookupOrDefault("buReleaseThresholdHBS", 80000);
        temperatureReleaseThresholdHBS_ = 
        fissionGasReleaseDict_.lookupOrDefault("temperatureReleaseThresholdHBS", 1000);
        damageReleaseThreshold_ = 
        fissionGasReleaseDict_.lookupOrDefault("damageReleaseThreshold", 0.85);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fgrSCIANTIXRIA::~fgrSCIANTIXRIA()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fgrSCIANTIXRIA::correct()
{    
    intragranularGasSwelling_.storePrevIter();
    intergranularGasSwelling_.storePrevIter();

    const scalarField& Ti = mesh_.lookupObject<volScalarField>("T").internalField();
    const scalarField& BuEffi = Effective_burn_up_;
    const scalarField& damage_ = mesh_.lookupObject<volScalarField>("damage").internalField();


    // Fuel cell volumes [m^3]
    const scalarField& V = mesh_.V();

    // Put the number of helium and gas moles released to zero
    scalar molFGRprev_(molFGR_);
    molHe_ *= 0;
    molFGR_ *= 0;

    const labelListList& matAddrList(mat_.matAddrList());

    // As SCIANTIX variables are scalarField and don't have oldTime
    // here we don't actually modify this variables.
    // They are going to be finally updated in the updateVariables function
    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                scalar gasReleased = 0.0;
                scalar HeReleased = 0.0;

                // Apply release criterion:

                // Release of the total gas inventory of the HBS based on a temperature criterion
                // The HBSThreshold_ is the value of BurnUp at which we consider that the HBS exists
                // The HBSReleaseThreshold_ is the considered temperature threshold  for release
                if 
                (
                    releaseHBS_ 
                    and (BuEffi[cellI]*1000 > buReleaseThresholdHBS_)
                    and (Ti[cellI] > temperatureReleaseThresholdHBS_)
                )
                {
                    gasReleased = 
                    Gas_released_[cellI] + Gas_boundary_[cellI] + Gas_grain_[cellI];
                    
                    intergranularGasSwelling_[cellI] *= 0.0;
                    intragranularGasSwelling_[cellI] *= 0.0;

                    HeReleased =
                    Helium_released_[cellI] + Helium_boundary_[cellI] + Helium_grain_[cellI];
                }
                else if (damage_[cellI] > damageReleaseThreshold_)
                {
                    gasReleased = 
                    Gas_released_[cellI] + Gas_boundary_[cellI];
                    
                    intergranularGasSwelling_[cellI] *= 0.0;

                    HeReleased =
                    Helium_released_[cellI] + Helium_boundary_[cellI];
                }
                else
                {   
                    gasReleased = Gas_released_[cellI];
                    intergranularGasSwelling_[cellI] = intergranularGasSwelling_.oldTime()[cellI];
                    intragranularGasSwelling_[cellI] = intragranularGasSwelling_.oldTime()[cellI];
                    HeReleased = Helium_released_[cellI];
                }

                // Relax swelling
                intragranularGasSwelling_[cellI] = 
                (
                    relax_*intragranularGasSwelling_[cellI] 
                    + (1- relax_)*intragranularGasSwelling_.prevIter()[cellI]
                );

                intergranularGasSwelling_[cellI] = 
                (
                    relax_*intergranularGasSwelling_[cellI] 
                    + (1- relax_)*intergranularGasSwelling_.prevIter()[cellI]
                );

                //***** Helium released *****//
                molHe_ += (HeReleased - Helium_released_[cellI])*V[cellI]/6.02e23;  

                //**** Total gas released *****//
                molFGR_ += (gasReleased - Gas_released[cellI])*V[cellI]/6.02e23;   
            }
        }
    }

    reduce(molFGR_, sumOp<scalar>());

    // Relax fgr
    molFGR_ = (1- relax_)*molFGRprev_ + relax_*molFGR_;
}

void Foam::fgrSCIANTIXRIA::updateVariables()
{
    const scalarField& Ti = mesh_.lookupObject<volScalarField>("T").internalField();
    const scalarField& BuEffi = Effective_burn_up_;
    const scalarField& damage_ = mesh_.lookupObject<volScalarField>("damage").internalField();

    // Fuel cell volumes [m^3]
    const scalarField& V = mesh_.V();

    // Put the number of helium and gas moles released to zero
    scalar molFGRprev_(molFGR_);
    molHe_ *= 0;
    molFGR_ *= 0;

    const labelListList& matAddrList(mat_.matAddrList());

    // As the SCIANTIX variables are scalarField and there is no oldTime
    // here we update the variables at the end of time step
    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                // Apply release criterion:

                // Release of the total gas inventory of the HBS based on a temperature criterion
                // The HBSThreshold_ is the value of BurnUp at which we consider that the HBS exists
                // The HBSReleaseThreshold_ is the considered temperature threshold  for release
                if 
                (
                    releaseHBS_ 
                    and (BuEffi[cellI]*1000 > buReleaseThresholdHBS_)
                    and (Ti[cellI] > temperatureReleaseThresholdHBS_)
                )
                {
                    Gas_released_[cellI] += Gas_boundary_[cellI] + Gas_grain_[cellI];
                    Gas_boundary_[cellI] = 0.0;
                    Gas_grain_[cellI] = 0.0;
                    intergranularGasSwelling_[cellI] *= 0.0;
                    intragranularGasSwelling_[cellI] *= 0.0;
                    Helium_released_[cellI] += Helium_boundary_[cellI] + Helium_grain_[cellI];
                    Helium_boundary_[cellI] = 0.0;
                    Helium_grain_[cellI] = 0.0;
                }
                else if (damage_[cellI] > damageReleaseThreshold_)
                {
                    Gas_released_[cellI] += Gas_boundary_[cellI];
                    Gas_boundary_[cellI] = 0.0;
                    intergranularGasSwelling_[cellI] *= 0.0;
                    Helium_released_[cellI] += Helium_boundary_[cellI];
                    Helium_boundary_[cellI] = 0.0;
                }
            }
        }
    }

    scalarField produced(mesh_.V()*pos(grainRadius_)*Gas_produced_);
    scalarField released(mesh_.V()*pos(grainRadius_)*Gas_released_);
    scalarField released0((mesh_.V()*Gas_released_t0_*pos(grainRadius_)));

    // Update fgr and fgr0 (for time stepping)
    fgr0_ = fgr_;
    fgr_ = (released - released0)/max(produced - released0, SMALL)*100;

    // scalar HeProduced(gSum(mesh_.V()*Helium_produced_*pos(grainRadius_)));
    // scalar HeReleased(gSum(mesh_.V()*Helium_released_*pos(grainRadius_)));

    // For debugging
    Info << "fission gas (Xe/Kr) atoms produced " << gSum(produced) << endl;
    Info << "fission gas (Xe/Kr) atoms released " << gSum(released) << endl;

    Info << "fgr% " << gSum(released)/max(gSum(produced),SMALL)*100 << endl;
    Info << "fgr% (RIA, relative) " << gSum(released - released0)/max(
        gSum(produced - released0),SMALL)*100 << nl << endl;
}


const Foam::wordList Foam::fgrSCIANTIXRIA::gasComponents() const
{
        wordList gasComp(3);
        gasComp[0] = "Xe";
        gasComp[1] = "Kr";
        gasComp[2] = "He";

        return gasComp;
}


const Foam::scalarList Foam::fgrSCIANTIXRIA::gasMols() const
{
        scalar XeMols_ = (0.268/0.301)*molFGR_;
        scalar KrMols_ = (0.033/0.301)*molFGR_;
        scalar HeMols_ = molHe_;
        
        // reduce(XeMols_, sumOp<scalar>());      
        // reduce(KrMols_, sumOp<scalar>());   
        
        scalarList gasMol(3, 0.0);
        
        gasMol[0] = XeMols_;
        gasMol[1] = KrMols_;
        gasMol[2] = HeMols_;
        
        return gasMol;
}


bool Foam::fgrSCIANTIXRIA::writeData(Ostream& os) const
{    
    const globalOptions& globalOpt
    (mesh_.lookupObject<globalOptions>("globalOptions"));

    // Obtain model angular fraction to calculate the release for the 
    // corresponding 360 degrees 3-D rod
    scalar angularFraction = globalOpt.angularFraction();

    // Calculate fgr%, print to screen and write in dictionary
    scalar produced(gSum(mesh_.V()*Gas_produced_*pos(grainRadius_)));
    scalar released(gSum(mesh_.V()*Gas_released_*pos(grainRadius_)));
    scalar HeProduced(gSum(mesh_.V()*Helium_produced_*pos(grainRadius_)));
    scalar HeReleased(gSum(mesh_.V()*Helium_released_*pos(grainRadius_)));

    // Write to dictionary
    os.writeKeyword("fgr produced (in the OFFBEAT model)") 
    << produced << token::END_STATEMENT << nl << nl;
    os.writeKeyword("fgr produced (in the corresponding 360 degree rod)") 
    << produced/angularFraction << token::END_STATEMENT << nl << nl;
    os.writeKeyword("fgr released (in the OFFBEAT model)") 
    << released << token::END_STATEMENT << nl << nl;
    os.writeKeyword("fgr released (in the corresponding 360 degree rod)") 
    << released/angularFraction << token::END_STATEMENT << nl << nl;

    os.writeKeyword("fgr%") 
    << released/max(produced,SMALL)*100 << token::END_STATEMENT << nl << nl;
      
    // Calculate and write relative fgr
    scalar deltaProduced(gSum(mesh_.V()*(
        Gas_produced_ -  Gas_released_t0_)*pos(grainRadius_)));
    scalar deltaReleased(gSum(mesh_.V()*(
        Gas_released_ - Gas_released_t0_)*pos(grainRadius_)));
    os.writeKeyword("fgr% (RIA, relative)") 
    << deltaReleased/max(deltaProduced,SMALL)*100 << token::END_STATEMENT << nl << nl;

    os.writeKeyword("He produced (in the OFFBEAT model)") 
    << HeProduced << token::END_STATEMENT << nl << nl;
    os.writeKeyword("He produced (in the corresponding 360 degree rod)") 
    << HeProduced/angularFraction << token::END_STATEMENT << nl << nl;
    os.writeKeyword("He released (in the OFFBEAT model)") 
    << HeReleased << token::END_STATEMENT << nl << nl;
    os.writeKeyword("He released (in the corresponding 360 degree rod)") 
    << HeReleased/angularFraction << token::END_STATEMENT << nl << nl;
    os.writeKeyword("He%") 
    << HeReleased/max(HeProduced,SMALL)*100 << token::END_STATEMENT << nl << nl;

    // If addToRegistry is true, avoid writing also in dictionary
    if(!addToRegistry_)
    {
        // os.writeKeyword("grainRadius_") 
        // << grainRadius_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_produced") 
        << Gas_produced_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_grain")
        << Gas_grain_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_grain_solution") 
        << Gas_grain_solution_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_grain_bubbles") 
        << Gas_grain_bubbles_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_boundary") 
        << Gas_boundary_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_released") 
        << Gas_released_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Gas_released_t0") 
        << Gas_released_t0_ << token::END_STATEMENT << nl << nl;

        os.writeKeyword("SCIANTIX::Effective_burn_up") 
        << Effective_burn_up_ << token::END_STATEMENT << nl << nl;
        
        os.writeKeyword("SCIANTIX::Helium_produced") 
        << Helium_produced_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Helium_grain") 
        << Helium_grain_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Helium_grain_solution") 
        << Helium_grain_solution_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Helium_grain_bubbles") 
        << Helium_grain_bubbles_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Helium_boundary") 
        << Helium_boundary_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Helium_released")         
        << Helium_released_ << token::END_STATEMENT << nl << nl;
    }
    
    return os.good();
}
// ************************************************************************* //
