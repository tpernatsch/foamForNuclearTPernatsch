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

#include "swellingParfumePyC.H"
#include "addToRunTimeSelectionTable.H"
#include "scalarFieldFieldINew.H"
#include "InterpolateTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingParfumePyC, 0);
    addToRunTimeSelectionTable(swellingModel, swellingParfumePyC, dictionary);

    const char* swellingParfumePyC::group_ = 
        ("behavioralModels::swelling::PyCParfume");

    defineParameter(swellingParfumePyC, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingParfumePyC, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingParfumePyC::swellingParfumePyC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    swellingModel(mesh, dict, baseDict),
    densityName_(dict.lookupOrDefault<word>("densityName", "rho" )),
    rho_(nullptr),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    fastFluxName_(dict.lookupOrDefault<word>("fastFluxName", "fastFlux")),
    fastFlux_(nullptr),
    BAF_(dict.lookupOrDefault("asFabricatedAnisotropy", 1.0)),
    sp_(dict.lookupOrDefault("sphereCoordinate", true)),
    fcf_(dict.lookupOrDefault("fluxConversionFactor", 1.0))
{
    //- Data file
    #include "swellingParfumePyCdata.H"

    //- Interpolation and out of bound method
    Method_ = interpolateTableBase::interpolationMethod::LINEAR;
    obMethod_ex_ = interpolateTableBase::outOfBoundsMethod::EXTRAPOLATION;
    obMethod_fix_ = interpolateTableBase::outOfBoundsMethod::FIXED;

    //- Read the temperature value
    tempValues_ = scalarField(dic.lookup("tempValues"));
    //- Read the BAF(BaconAnisotropyFactor) value for radial and tangential components
    BAFValues_r_ = scalarField(dic.lookup("BAFValues_r"));
    BAFValues_t_ = scalarField(dic.lookup("BAFValues_t"));
    //- Read density value for 1D table of strain_iso and strain_r - strain_t (strain_r_t)
    densityValue_ = scalarField(dic.lookup("densityValues"));
    strain_iso_ = scalarField(dic.lookup("strain_iso"));
    strain_r_t_ = scalarField(dic.lookup("strain_r_t"));


    //- Read the Polynomial coefficients data
    //- Raidal components
    a1_r_ = PtrList<scalarField>(dic.lookup("a1_r"), scalarFieldFieldINew());
    a2_r_ = PtrList<scalarField>(dic.lookup("a2_r"), scalarFieldFieldINew());
    a3_r_ = PtrList<scalarField>(dic.lookup("a3_r"), scalarFieldFieldINew());
    a4_r_ = PtrList<scalarField>(dic.lookup("a4_r"), scalarFieldFieldINew());
    //- Tangential components
    a1_t_ = PtrList<scalarField>(dic.lookup("a1_t"), scalarFieldFieldINew());
    a2_t_ = PtrList<scalarField>(dic.lookup("a2_t"), scalarFieldFieldINew());
    a3_t_ = PtrList<scalarField>(dic.lookup("a3_t"), scalarFieldFieldINew());
    a4_t_ = PtrList<scalarField>(dic.lookup("a4_t"), scalarFieldFieldINew());

    //- Correspond the Polynomial coefficients data with temperature.Values at
    //- temperatures outside validity range are extrapolated from values at 600
    //- and 1032(below range) or 1032 and 1350 (above range), respectively.

    //- Radial components
    table_a1_r_.set(new scalarFieldInterpolateTable(tempValues_, a1_r_, Method_, obMethod_ex_));
    table_a2_r_.set(new scalarFieldInterpolateTable(tempValues_, a2_r_, Method_, obMethod_ex_));
    table_a3_r_.set(new scalarFieldInterpolateTable(tempValues_, a3_r_, Method_, obMethod_ex_));
    table_a4_r_.set(new scalarFieldInterpolateTable(tempValues_, a4_r_, Method_, obMethod_ex_));
    //- Tangential components
    table_a1_t_.set(new scalarFieldInterpolateTable(tempValues_, a1_t_, Method_, obMethod_ex_));
    table_a2_t_.set(new scalarFieldInterpolateTable(tempValues_, a2_t_, Method_, obMethod_ex_));
    table_a3_t_.set(new scalarFieldInterpolateTable(tempValues_, a3_t_, Method_, obMethod_ex_));
    table_a4_t_.set(new scalarFieldInterpolateTable(tempValues_, a4_t_, Method_, obMethod_ex_));
    //- Taken table_a1_r as expample:
    /*
    600  | (-1.24080 -1.10640 -0.94333 -0.78045 -0.15714  0.40265)
    1032 | (-1.52390 -2.07520 -2.00470 -1.81690 -1.18540 -0.45900)
    1350 | (-1.42840 -1.54330 -1.49640 -0.89522  1.20930  3.71620)
    */

    //- 1D table for epsilon.
    //- Swelling strain at current density of 1100 C and fast fluence of 3.7x1025 n/m2 throng the Interpolation
    //- Values at densities below 1.0 g/cm3 are extrapolated from values at 1.2 and 1.0 g/cm33
    //- Values at densities above 2.0 g/cm33 are extrapolated from values at 1.96 and 2.0 g/cm3.
    table_epsilon_iso_.set(new scalarInterpolateTable(densityValue_, strain_iso_, Method_, obMethod_ex_));
    table_epsilon_r_t_.set(new scalarInterpolateTable(densityValue_, strain_r_t_, Method_, obMethod_ex_));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingParfumePyC::~swellingParfumePyC()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingParfumePyC::correct
(
    const scalarField& T,
    const labelList& addr
)
{
    // Check which F and delta should be used
    const dimensionedScalar& F =
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() :
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta =
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() :
    swellingModel::delta_epsilonSwelling();

    const fvMesh& mesh_ = epsilonSwelling_.mesh();

    if(fastFluence_ == nullptr or rho_==nullptr or fastFlux_ == nullptr)
    {
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
        fastFlux_ = &mesh_.lookupObject<volScalarField>(fastFluxName_);
        rho_ = &mesh_.lookupObject<volScalarField>(densityName_);
    }

    const symmTensorField& swellingIOld = epsilonSwelling_.oldTime().internalField();

    symmTensorField& swellingI = epsilonSwelling_.ref();

    const volScalarField& phii = *fastFluence_;
    const volScalarField& rhooo = *rho_;
    const volScalarField& fluxx = *fastFlux_;
    //Only read the inital density
    static scalarField rhoo = rhooo;

    //- Reference to the bulit table
    const scalarFieldInterpolateTable& table_a1_r_ref(table_a1_r_());
    const scalarFieldInterpolateTable& table_a2_r_ref(table_a2_r_());
    const scalarFieldInterpolateTable& table_a3_r_ref(table_a3_r_());
    const scalarFieldInterpolateTable& table_a4_r_ref(table_a4_r_());

    const scalarFieldInterpolateTable& table_a1_t_ref(table_a1_t_());
    const scalarFieldInterpolateTable& table_a2_t_ref(table_a2_t_());
    const scalarFieldInterpolateTable& table_a3_t_ref(table_a3_t_());
    const scalarFieldInterpolateTable& table_a4_t_ref(table_a4_t_());

    const scalarInterpolateTable& table_epsilon_iso_ref(table_epsilon_iso_());
    const scalarInterpolateTable& table_epsilon_r_t_ref(table_epsilon_r_t_());

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        //- The model requires 1e25 n/m^2, it is bounded at 3.96e25 n/m^2
        const scalar phi = min(phii[cellI]*1e4/1e25*fcf_, 3.96);
        //- Fast flux in unit 1e25n/m²/s
        const scalar flux = fluxx[cellI]*1e4/1e25*fcf_;
        //- Density in g/cm3
        const scalar rho = rhoo[cellI]/1e3;
        //- Temperature in Celsius
        const scalar Ti = T[cellI] - 273.15;
        //- Do the interpolation for the temperature
        scalarField a1_r_T(table_a1_r_ref(Ti));
        scalarField a2_r_T(table_a2_r_ref(Ti));
        scalarField a3_r_T(table_a3_r_ref(Ti));
        scalarField a4_r_T(table_a4_r_ref(Ti));
        scalarField a1_t_T(table_a1_t_ref(Ti));
        scalarField a2_t_T(table_a2_t_ref(Ti));
        scalarField a3_t_T(table_a3_t_ref(Ti));
        scalarField a4_t_T(table_a4_t_ref(Ti));

        //- Correspond the scalarField above with the BAF data.
        //- BAFvalue outside validity range are clamped at the bounds.
        scalarInterpolateTable a1_r_T_BAF(BAFValues_r_, a1_r_T, Method_, obMethod_fix_);
        scalarInterpolateTable a2_r_T_BAF(BAFValues_r_, a2_r_T, Method_, obMethod_fix_);
        scalarInterpolateTable a3_r_T_BAF(BAFValues_r_, a3_r_T, Method_, obMethod_fix_);
        scalarInterpolateTable a4_r_T_BAF(BAFValues_r_, a4_r_T, Method_, obMethod_fix_);
        scalarInterpolateTable a1_t_T_BAF(BAFValues_t_, a1_t_T, Method_, obMethod_fix_);
        scalarInterpolateTable a2_t_T_BAF(BAFValues_t_, a2_t_T, Method_, obMethod_fix_);
        scalarInterpolateTable a3_t_T_BAF(BAFValues_t_, a3_t_T, Method_, obMethod_fix_);
        scalarInterpolateTable a4_t_T_BAF(BAFValues_t_, a4_t_T, Method_, obMethod_fix_);

        //- Parameters
        scalar a1_iso, a2_iso, a3_iso, a4_iso;
        scalar a1_r, a2_r, a3_r, a4_r;
        scalar a1_t, a2_t, a3_t, a4_t;
        //- Interpolation for "a"
        a1_iso = a1_r_T_BAF(1.0);
        a2_iso = a2_r_T_BAF(1.0);
        a3_iso = a3_r_T_BAF(1.0);
        a4_iso = a4_r_T_BAF(1.0);
        a1_r = a1_r_T_BAF(BAF_);
        a2_r = a2_r_T_BAF(BAF_);
        a3_r = a3_r_T_BAF(BAF_);
        a4_r = a4_r_T_BAF(BAF_);
        a1_t = a1_t_T_BAF(BAF_);
        a2_t = a2_t_T_BAF(BAF_);
        a3_t = a3_t_T_BAF(BAF_);
        a4_t = a4_t_T_BAF(BAF_);

        //- epsilon value at density 1.96 g/cm3, the strain is in %.
        // const scalar epsilon_iso_rho0 = a1_iso*phi+a2_iso*pow(phi,2)+a3_iso*pow(phi,3)+a4_iso*pow(phi,4);
        // const scalar epsilon_r_rho0 = a1_r*phi+a2_r*pow(phi,2)+a3_r*pow(phi,3)+a4_r*pow(phi,4);
        // const scalar epsilon_t_rho0 = a1_t*phi+a2_t*pow(phi,2)+a3_t*pow(phi,3)+a4_t*pow(phi,4);
        const scalar epsilon_iso_rho0 = a1_iso+2.0*a2_iso*phi+3.0*a3_iso*pow(phi,2)+4.0*a4_iso*pow(phi,3);
        const scalar epsilon_r_rho0 = a1_r+2.0*a2_r*phi+3.0*a3_r*pow(phi,2)+4.0*a4_r*pow(phi,3);
        const scalar epsilon_t_rho0 = a1_t+2.0*a2_t*phi+3.0*a3_t*pow(phi,2)+4.0*a4_t*pow(phi,3);

        //- Isotropic scaling factor μ.
        const scalar mu = table_epsilon_iso_ref(rho) / table_epsilon_iso_ref(1.96);
        //- Anisotropic scaling factor ν.
        const scalar nu = table_epsilon_r_t_ref(rho) / table_epsilon_r_t_ref(1.96);

        //- epsilon value at current density
        const scalar epsilon_r = (mu*epsilon_iso_rho0 + nu*(epsilon_r_rho0-epsilon_iso_rho0))/1e2;
        const scalar epsilon_t = (mu*epsilon_iso_rho0 + nu*(epsilon_t_rho0-epsilon_iso_rho0))/1e2;

        //- Strain in spherical coordinate. Transverse strain is firstly given
        symmTensor nominalValue_sphere = (epsilon_t*flux*mesh_.time().deltaT().value())* I;
        //- Raidal strain replace the value of transverse strain at yy and zz.
        nominalValue_sphere.xx() = epsilon_r*flux*mesh_.time().deltaT().value();

        // nominalValue_sphere.yy() = epsilon_t*flux*mesh_.time().deltaT().value();
        // nominalValue_sphere.zz() = epsilon_t*flux*mesh_.time().deltaT().value();

        if (sp_ == true)
        {
          //Perturb for sensitivity analysis
          nominalValue_sphere = nominalValue_sphere * F.value() 
                            + delta.value() * symmTensor::one; 
          
          swellingI[cellI] = swellingIOld[cellI]+nominalValue_sphere;
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

          //- Calculate the angels theta and phi in sphere coordinate
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
          symmTensor nominalValue = symm(T1 & nominalValue_sphere & T2);
          //*************************coordinate conversion end***********************
          
          //Perturb for sensitivity analysis
          nominalValue = nominalValue * F.value() + delta.value() * symmTensor::one;

          swellingI[cellI] = swellingIOld[cellI] + nominalValue;
        }
    }
}


// ************************************************************************* //
