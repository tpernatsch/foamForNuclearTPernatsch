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

#include "swellingCorrelationPyC.H"
#include "addToRunTimeSelectionTable.H"

#include "scalarFieldFieldINew.H"
#include "InterpolateTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingCorrelationPyC, 0);
    addToRunTimeSelectionTable(swellingModel, swellingCorrelationPyC, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingCorrelationPyC::swellingCorrelationPyC
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    swellingModel(mesh, dict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    sp_(dict.lookupOrDefault("sphereCoordinate", true)),
    fcf_(dict.lookupOrDefault("fluxConversionFactor", 1.0)),
    Ar_(dict.lookup("radialCoefficients")),
    At_(dict.lookup("tangentialCoefficients")),
    perturb(1.0)
{
    if(dict.found("swelling"))
    {
        const dictionary& swellingDict = dict.subDict("swelling");
        perturb = swellingDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingCorrelationPyC::~swellingCorrelationPyC()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingCorrelationPyC::correct
(
    const scalarField& T,
    const labelList& addr
)
{
    const fvMesh& mesh_ = epsilonSwelling_.mesh();

    if(fastFluence_ == nullptr)
    {
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();

    const volScalarField& phii = *fastFluence_;

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        //- The model requires 1e25 n/m^2,
        const scalar phi = phii[cellI]*1e4/1e25*fcf_;

        const scalar epsilon_r = Ar_[0]*phi + Ar_[1]/2.0*phi*phi + Ar_[2]/3.0*pow(phi,3)
                               + Ar_[3]/4.0*pow(phi,4) + Ar_[4]/5.0*pow(phi,5)
                               + Ar_[5]/6.0*pow(phi,6);
        const scalar epsilon_t = At_[0]*phi + At_[1]/2.0*phi*phi + At_[2]/3.0*pow(phi,3)
                               + At_[3]/4.0*pow(phi,4) + At_[4]/5.0*pow(phi,5)
                               + At_[5]/6.0*pow(phi,6);

        //- Strain in spherical coordinate. Radial strain is firstly given
        symmTensor nominalValue_sphere = epsilon_r * I;
        //- Transverse strain replace the value of raidal strain at yy and zz.
        nominalValue_sphere.yy() = epsilon_t;
        nominalValue_sphere.zz() = epsilon_t;

        if (sp_ == true)
        {
          swellingI[cellI] = nominalValue_sphere*perturb;
        }
        else
        {
          //*************************coordinate conversion***********************
          //- obtain the coordinate of elements
          const scalar x_coord = mesh_.C()[cellI].x();
          const scalar y_coord = mesh_.C()[cellI].y();
          const scalar z_coord = mesh_.C()[cellI].z();
          //- Radius
          const scalar R = sqrt(pow(x_coord,2)+pow(y_coord,2)+pow(z_coord,2));

          scalar theta, phi_angle;

          //- Calculate the theta and phi in sphere coordinate
          if (R <= 1e-10) //- For the element where its center is the center of the sphere.
          {
            theta = 0.0;
            phi_angle = 0.0;
          }
          else
          {
            theta = acos(z_coord/R);
            phi_angle = atan(y_coord/x_coord);
          }
          //- Transformation matrix. According to :
          //- https://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
          const tensor T1(sin(theta)*cos(phi_angle), cos(theta)*cos(phi_angle), -sin(phi_angle),
                          sin(theta)*sin(phi_angle), cos(theta)*sin(phi_angle), cos(phi_angle),
                          cos(theta),          -sin(theta),         0.0);
          const tensor T2(sin(theta)*cos(phi_angle), sin(theta)*sin(phi_angle), cos(theta),
                          cos(theta)*cos(phi_angle), cos(theta)*sin(phi_angle), -sin(theta),
                          -sin(phi_angle),           cos(phi_angle),            0.0);

          //inner product
          const symmTensor nominalValue = symm(T1 & nominalValue_sphere & T2);
          //*************************coordinate conversion end***********************

          swellingI[cellI] = nominalValue*perturb;
        }


    }

    // Check which F and delta should be used
    const dimensionedScalar& F =
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() :
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta =
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() :
    swellingModel::delta_epsilonSwelling();

    // Perturb epsilon swelling for sensitivity analysis
    applyParameters(epsilonSwelling_.ref(), addr, F, delta);

}


// ************************************************************************* //
