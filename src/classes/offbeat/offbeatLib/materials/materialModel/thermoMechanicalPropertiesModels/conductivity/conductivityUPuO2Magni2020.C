
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

#include "conductivityUPuO2Magni2020.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2Magni2020, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2Magni2020, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Magni2020::conductivityUPuO2Magni2020
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    kappaInf_(1.755),     
    phi_(128.75e3),     
    OM_(dict.lookupOrDefault<scalar>("oxygenMetalRatio",2)),
    x_(2-OM_),     
    Pu_(),
    ratioUMetal_(readScalar(dict.lookup("ratioUMetal"))),
    ratioPuMetal_(readScalar(dict.lookup("ratioPuMetal"))),
#ifdef OPENFOAMFOUNDATION
    isotopesU_
    (
        dict.lookupOrDefault<scalarList>
        (
            "isotopesU", scalarList(4,0.0)
        )
    ),
    isotopesPu_
    (
        dict.lookupOrDefault<scalarList>
        (
            "isotopesPu", scalarList(5,0.0)
        )
    ),
#elif OPENFOAMESI
    isotopesU_
    (
        dict.getOrDefault<scalarList>
        (
            "isotopesU", scalarList(4,0.0)
        )
    ),
    isotopesPu_
    (
        dict.getOrDefault<scalarList>
        (
            "isotopesPu", scalarList(5,0.0)
        )
    ),
#endif
    porosity_(1-readScalar(dict.lookup("densityFraction"))),
    A0_(0.01926),     
    Ax_(1.06e-6),     
    APu_(2.63e-8),     
    B0_(2.39e-4),     
    BPu_(1.37e-13),     
    D_(5.27e9),     
    E_(17109.5),
    perturb(1)  
{
    if(dict.found("conductivity"))
    {
        const dictionary& conductivityDict = dict.subDict("conductivity");

        kappaInf_ = conductivityDict.lookupOrDefault<scalar>("kappaInf", 1.755);
        phi_ = conductivityDict.lookupOrDefault<scalar>("phi", 128.75e3);
        A0_ = conductivityDict.lookupOrDefault<scalar>("A0", 0.01926);
        Ax_ = conductivityDict.lookupOrDefault<scalar>("Ax", 1.06e-6);
        APu_ = conductivityDict.lookupOrDefault<scalar>("APu", 2.63e-8);
        B0_ = conductivityDict.lookupOrDefault<scalar>("B0", 2.39e-4);
        BPu_ = conductivityDict.lookupOrDefault<scalar>("BPu", 1.37e-13);
        D_ = conductivityDict.lookupOrDefault<scalar>("D", 5.27e9);
        E_ = conductivityDict.lookupOrDefault<scalar>("E", 17109.5);

        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

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

    Info << "\tMolar mass U is : "  << MU  << " g/mol" << endl;
    Info << "\tMolar mass Pu is : " << MPu << " g/mol" << endl;

    // Compute atomic concentration of Pu in HM:
    Pu_ = ratioPuMetal_/MPu * pow(ratioPuMetal_/MPu+ratioUMetal_/MU, -1);

    Info << "\tAtomic concentation of Pu is : " << Pu_ << " at.%" << endl;

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Magni2020::~conductivityUPuO2Magni2020()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2Magni2020::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }
 
    const scalarField& Bu = Bu_->internalField();

    //- Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];              

        const scalar Bui = Bu[cellI];

        const scalar Ti = T[cellI];

        //- Computing k0 (= thermal conductivity at zero Bu )            
        double term1 = 1/(A0_ + Ax_*x_ + APu_*Pu_ + (B0_ + BPu_*Pu_)*Ti);
        double term2 = D_/pow(Ti,2.0)*exp(-E_/Ti);
        double term3 = pow(1-porosity_ , 2.5);
        double k0    = (term1 + term2)*term3;

        //- Computing kirr (= thermal conductivity as a f(Bu) )
        const scalar nominalValue = kappaInf_ + (k0 - kappaInf_)*exp(-Bui/phi_);

        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
