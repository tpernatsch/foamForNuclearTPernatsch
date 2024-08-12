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

#include "thermalExpansionPARFUMEPyC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionPARFUMEPyC, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel,
        thermalExpansionPARFUMEPyC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionPARFUMEPyC::thermalExpansionPARFUMEPyC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    BAF_(dict.lookupOrDefault("asFabricatedAnisotropy", 1.0)),
    sp_(dict.lookupOrDefault("sphereCoordinate", true)),
    par1(30.0),
    par2(37.5),
    par3(0.11),
    par4(673.0),
    par5(700.0),
    par6(36.0),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", 30.0);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 37.5);
        par3 = thermalExpansionDict.lookupOrDefault<scalar>("par3", 0.11);
        par4 = thermalExpansionDict.lookupOrDefault<scalar>("par4", 673.0);
        par5 = thermalExpansionDict.lookupOrDefault<scalar>("par5", 700.0);
        par6 = thermalExpansionDict.lookupOrDefault<scalar>("par6", 36.0);

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionPARFUMEPyC::~thermalExpansionPARFUMEPyC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionPARFUMEPyC::correct
(
    symmTensorField& sf,
    const scalarField& T,
    const labelList& addr
)
const
{
    //- Orientation parameters
    const scalar Rr = 2.0 / (2.0 + BAF_);
    const scalar Rt = (1.0 + BAF_) / (2.0 + BAF_);

    forAll(addr, i)
    {
        const label cellI = addr[i];

        // Temperature must be given in K
        const scalar Ti = T[cellI];
        const scalar Tref = Tref_.value();
        //- radial thermal expansion coef
        const scalar alpha_r = (par1-par2*Rr)*(1.0+par3*(Ti-par4)/par5)*1e-6;
        //- tagential thermal expansion coef
        const scalar alpha_t = (par6*pow(Rt-1.0,2.0)+1.0)*(1.0+par3*(Ti-par4)/par5)*1e-6;

        //-Thermal strain in spherical coordinate. Radial thermal strain is firstly given
        symmTensor nominalValue_sphere = alpha_r*(Ti-Tref)*I;
        symmTensor nominalValue;

        //- Transverse thermal strain replace the value of raidal strain at yy and zz.
        nominalValue_sphere.yy() = alpha_t*(Ti-Tref);
        nominalValue_sphere.zz() = alpha_t*(Ti-Tref);

        // NOTE: avoids instabilities when trying to simulate a material at
        // constant temperature equal to Tref
        if
        (
          nominalValue_sphere.xx() < 1e-7 &&
          nominalValue_sphere.yy() < 1e-7 &&
          nominalValue_sphere.zz() < 1e-7
        )
        {
            nominalValue_sphere *= 0;
        }

        if (sp_ == true)
        {
          sf[cellI] = nominalValue_sphere*perturb;
        }
        else //in case of non-1D case
        {
          //*************************coordinate conversion***********************
          //- obtain the coordinate of elements
          const scalar x_coord = mesh_.C()[cellI].x();
          const scalar y_coord = mesh_.C()[cellI].y();
          const scalar z_coord = mesh_.C()[cellI].z();
          //- Radius
          const scalar R = sqrt(pow(x_coord,2)+pow(y_coord,2)+pow(z_coord,2));

          //- Calculate the angels theta and phi in sphere coordinate
          const scalar theta = acos(z_coord/R);
          const scalar phi = atan(y_coord/x_coord);
          //- Transformation matrix. According to :
          //- https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
          const tensor T1(sin(theta)*cos(phi), cos(theta)*cos(phi), -sin(phi),
                          sin(theta)*sin(phi), cos(theta)*sin(phi), cos(phi),
                          cos(theta),          -sin(theta),         0.0);
          const tensor T2(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta),
                          cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta),
                          -sin(phi),           cos(phi),            0.0);

          //inner product
          nominalValue = symm(T1 & nominalValue_sphere & T2);

          sf[cellI] = nominalValue*perturb;
        }

    }
}

// ************************************************************************* //
