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

#include "UPuO2meltingMagni2020.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPuO2meltingMagni2020, 0);
    addToRunTimeSelectionTable(failureModel, UPuO2meltingMagni2020, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPuO2meltingMagni2020::UPuO2meltingMagni2020
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    failureModel(mesh, dict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    TmUO2_(failureModelDict_.lookupOrDefault<scalar>("Tmelt", 3147.0)),
    Pu_(),
    ratioUMetal_(readScalar(dict.lookup("ratioUMetal"))),
    ratioPuMetal_(readScalar(dict.lookup("ratioPuMetal"))),
#ifdef OPENFOAMFOUNDATION
    isotopesU_(dict.lookup<scalarList>("isotopesU")),
    isotopesPu_(dict.lookup<scalarList>("isotopesPu")),
#elif OPENFOAMESI
    isotopesU_(dict.get<scalarList>("isotopesU")),
    isotopesPu_(dict.get<scalarList>("isotopesPu")),
#endif
    OM_(readScalar(dict.lookup("oxygenMetalRatio"))),
    x_(2-OM_),
    gammaPu_(364.85),
    gammaX_(1014.15),
    TmInf_(2964.92),
    delta_(40.43e3)
{

    // Molar mass U
    const scalar MU = 234 * isotopesU_[0]
                    + 235 * isotopesU_[1]
                    + 236 * isotopesU_[2]
                    + 238 * isotopesU_[3];

    // Molar mass Pu
    const scalar MPu = 238 * isotopesPu_[0]
                     + 239 * isotopesPu_[1]
                     + 240 * isotopesPu_[2]
                     + 241 * isotopesPu_[3]
                     + 242 * isotopesPu_[4];

    //Info << "\tMolar mass U is : "  << MU  << " g/mol" << endl;
    //Info << "\tMolar mass Pu is : " << MPu << " g/mol" << endl;

    // Compute atomic concentration of Pu in HM:
    Pu_ = ratioPuMetal_/MPu * pow(ratioPuMetal_/MPu+ratioUMetal_/MU, -1);

    //Info << "\tAtomic concentation of Pu is : " << Pu_ << " at.%" << endl;    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UPuO2meltingMagni2020::~UPuO2meltingMagni2020()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::UPuO2meltingMagni2020::checkFailure
(
    const labelList& addr
)
{
    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }
 
    const scalarField& Bu = Bu_->internalField();

    // Set failure switch to false
    failed_ = false;

    // Set counter of failed cells to 0
    scalar counter(0);

    // Compute Tm,0
    const scalar Tm0 = TmUO2_ - gammaPu_*Pu_ - gammaX_*x_;
    
    // Constant reference to the T field
    const volScalarField& T(mesh_.lookupObject<volScalarField>("T"));

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        // Compute Tm,irr
        const scalar TmIrr 
            = TmInf_ + (Tm0 - TmInf_)*Foam::exp(-Bu[cellI]/delta_);

        if ( T[cellI] >= TmIrr)
        {
            if( !mesh_.foundObject<volScalarField>("failedMaterial"))
            {
                //- Initialize failedMaterial field (from parent class "failureModel") 
                initializeFailedMaterialField();
            }

            failedMaterial_()[cellI] = 1;
            failed_ = true;
            counter++;
        }
    }

    if ( failed_ )
    {

        if ( stopIfFailed_ )
        {
            // Write fields files in the time-step directory at which the failure occurs
            mesh_.time().write();

            // Fatal error
            FatalErrorIn("Foam::UPuO2meltingMagni2020::checkFailure(const labelList& addr)") << nl
            << "    Failure occurred in " << counter << " out of " << addr.size() 
            << " cells of the material."  << nl      << exit(FatalError); 
        }

        //- Warning
        WarningIn("Foam::UPuO2meltingMagni2020::checkFailure(const labelList& addr)") << nl  
        << "    Failure occurred in " << counter << " out of " << addr.size() 
        << " cells of the material."  << nl      << endl;
    }
}


// ************************************************************************* //

