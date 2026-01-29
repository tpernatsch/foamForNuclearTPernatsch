/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "NaCoolantChannelfvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private functions * * * * * * * * * * * * * * * * //

Foam::scalarField NaCoolantChannelfvPatchScalarField::computeSliceHeights()
{
    // Get n slices of the cellZone to which this patch belong. The heights of 
    // the slices are used to get the dy to integrate the fluid enthalpy axially

    // Notice : this BC must then be used for a patch that belongs to a single
    //          cellzone (i.e. the slices of the cellZone must axially 'cover'
    //          the whole patch)

    // Ref. to mesh
    const fvMesh& mesh(patch().boundaryMesh().mesh());

    // Get name of cellzone to which this patch is attached
    const labelList faceCells(patch().faceCells());
    
    // Get cell id of one of the cells attached to this patch (e.g. the first)
    const label cellIdx(faceCells[0]);

    // Loop to find cellZone name
    word cellZoneName("");
    label cellZoneID(-1);
    forAll(mesh.cellZones(), i)
    {
        const labelList addr(mesh.cellZones()[i]);
        forAll(addr, j)
        {
            if( addr[j] == cellIdx)
            {
                cellZoneName = mesh.cellZones()[i].name();
                cellZoneID = i;
                break;
            }
        }
    }
    if( cellZoneName == "" )
    {
        FatalErrorIn("NaCoolantChannelfvPatchScalarField::computeSliceHeights()")
        << "Cannot find cellZone name for the patch " << patch().name()
        << nl << nl << exit(FatalError);
    }

    // Get list of all the addressings corresponding to ALL slices
    const labelListList allSlicesAddrList(mapper_->sliceAddrList());

    // Loop over all the slice addressings and keep only the ones belonging to
    // this cellzone
    labelListList sliceAddrList(0);
    forAll(allSlicesAddrList, i)
    {
        // Take the label of one of the cells of this slice (e.g. the first)
        const label cellI(allSlicesAddrList[i][0]);

        // Loop over addr of this cellZone to find out if cellI is in there
        forAll(mesh.cellZones()[cellZoneID], j)
        {
            if( cellI == mesh.cellZones()[cellZoneID][j])
            {
                sliceAddrList.append(allSlicesAddrList[i]);
            }
        }
    }

    const vectorField& points = patch().boundaryMesh().mesh().points();
    const labelListList& cellPoints = patch().boundaryMesh().mesh().cellPoints();
    vector pinDirection(0, 0, 1);
    scalarField z = (points & pinDirection);

    scalar nSlices(sliceAddrList.size());
    scalarField sliceHeights(nSlices, 0.0);

    forAll(sliceAddrList, sliceID)
    {
        const labelList addr(sliceAddrList[sliceID]);

        scalar zMax = -GREAT;
        scalar zMin = GREAT;
        
        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];
            const labelList& ptIds = cellPoints[cellI];
            
            forAll(ptIds, i)
            {
                const label ptI = ptIds[i];
                zMin = min(zMin, z[ptI]);
                zMax = max(zMax, z[ptI]);
            }
        }

        reduce(zMin, minOp<scalar>());
        reduce(zMax, maxOp<scalar>());

        sliceHeights[sliceID] = zMax - zMin;
    }

    return sliceHeights;
}

Foam::scalar NaCoolantChannelfvPatchScalarField::NaEnthalpy
(
    scalar T
)
{
    // Correlation from ANL/RE-92/2 report (Fink and Leibowitz, 1995)
    // Enthalpy relative to T = 298.15
    return -365770 + 1658.2*T - 0.42395*pow(T,2.0) + 1.4847e-4*pow(T,3.0) + 2992600/T;
}

Foam::scalar NaCoolantChannelfvPatchScalarField::NaConductivity
(
    scalar T
)
{
    if( modelKappa_ == "constant" )
    {
        // return 72;
        return constKappa_;
    }
    else if( modelKappa_ == "FinkLeibowitz" )
    {
        // Correlation from ANL/RE-92/2 report (Fink and Leibowitz, 1995)
        return 124.67 - 0.11381*T + 5.5226e-5*pow(T,2.0) - 1.1842e-8*pow(T,3.0);
    }
    else
    {
        FatalErrorInFunction()  
        << "Incorrectly specified Cp model name." << nl
        << "Available models are :" << nl  
        << "  - constant" << nl  
        << "  - FinkLeibowitz" << nl  
        << abort(FatalError);
    }
}

Foam::scalar NaCoolantChannelfvPatchScalarField::NaCp
(
    scalar T
)
{
    if( modelCP_ == "constant" )
    {
        // Constant Cp
        // return 1275;
        return constCp_;
    }
    else if( modelCP_ == "Vargaftic" )
    {
        // Vargaftic,1975 for Pure Na
        // Temperature Range 373-1500 K
        return 1732 - 1.266*T + 1.031e-3*pow(T,2.0) - 2.369e-7*pow(T,3.0);
    }
    else
    {
        FatalErrorInFunction()  
        << "Incorrectly specified Cp model name." << nl
        << "Available models are :" << nl  
        << "  - constant" << nl  
        << "  - Vargaftic" << nl  
        << abort(FatalError);
    }
}

Foam::scalar NaCoolantChannelfvPatchScalarField::NaVisc
(
    scalar T
)
{
    if( modelViscosity_ == "constant" )
    {
        // return 0.25e-3;
        return constVisc_;
    }
    else if( modelViscosity_ == "Shpilrain" )
    {
        // Shpil'rain, 1985 for Pure Na
        // Temperature Range 371-2400 K
        return exp( 556.835/T - 0.3958*log(T) - 6.4406 );
    }
    else
    {
        FatalErrorInFunction()  
        << "Incorrectly specified Cp model name." << nl
        << "Available models are :" << nl  
        << "  - constant" << nl  
        << "  - Shpilrain" << nl  
        << abort(FatalError);
    }
}

Foam::scalar NaCoolantChannelfvPatchScalarField::NaPrandtl
(
    scalar T
)
{
    // Compute Pr = cp*mu/k
    return ( NaCp(T) * NaVisc(T) ) / NaConductivity(T);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NaCoolantChannelfvPatchScalarField::
NaCoolantChannelfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    kappaName_("k"),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    pOverD_(0.0),
    A_(0.0),
    relaxHTC_(0.15),
    DH_(0.0),
    massFlowRate_(0.0),
    massFlowRateList_(),
    inletTemperature_(0.0),
    inletTemperatureList_(),
    modelKappa_("FinkLeibowitz"),
    constKappa_(0),
    modelCP_("Vargaftic"),
    constCp_(0),
    modelViscosity_("Shpilrain"),
    constVisc_(0),
    modelNusselt_("Schad"),
    constNusselt_(0),
    writeCoolantT_(false),
    Tcoolant_(nullptr),
    Hcoolant_(nullptr),
    heatFlux_(p.size(), 0.0),
    mapper_(nullptr)
{}

NaCoolantChannelfvPatchScalarField::
NaCoolantChannelfvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict,
     const bool valueRequired
 )
 :
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    kappaName_
    (
        dict.lookupOrDefault<word>("kappa", "k")
    ),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    pOverD_(0.0),
    A_(readScalar(dict.lookup("flowArea"))),
    relaxHTC_(dict.lookupOrDefault("relaxHTC", 0.15)),
    DH_(readScalar(dict.lookup("hydraulicDiameter"))),
    massFlowRate_(0.0),
    massFlowRateList_(),
    inletTemperature_(0.0),
    inletTemperatureList_(),
    modelKappa_(dict.lookupOrDefault<word>("conductivityModel","FinkLeibowitz")),
    constKappa_(0),
    modelCP_(dict.lookupOrDefault<word>("heatCapacityModel","Vargaftic")),
    constCp_(0),
    modelViscosity_(dict.lookupOrDefault<word>("viscosityModel","Shpilrain")),
    constVisc_(0),
    modelNusselt_(dict.lookupOrDefault<word>("NusseltModel","Lyon")),
    constNusselt_(0),
    writeCoolantT_(dict.lookupOrDefault<bool>("writeCoolantT",false)),
    Tcoolant_(nullptr),
    Hcoolant_(nullptr),
    heatFlux_(p.size(), 0.0),
    mapper_(nullptr)
{
    //- Read fluid mass flow rate (either list or fixed)
    if(dict.found("massFlowRateList"))
    {
#ifdef OPENFOAMFOUNDATION            
        massFlowRateList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "massFlowRateList", dict.subDict("massFlowRateList")
            )
        );
#elif OPENFOAMESI            
        massFlowRateList_.reset
        ( 
            Function1<scalar>::New
            (
                "massFlowRateList", dict.subDict("massFlowRateList")
            )
        );
#endif        
        
        massFlowRate_ = 
        massFlowRateList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
    else if (dict.found("massFlowRate"))
    {
        massFlowRate_ = readScalar(dict.lookup("massFlowRate"));
    }
    else
    {
        FatalErrorInFunction() 
            << "Coolant channel BC for " << patch().name() << " requires:" << nl
            << "- \"massFlowRateList\" or" << nl
            << "- \"massFlowRate\"" << abort(FatalError) << endl;
    }

        //- Read fluid inlet temperature (either list or fixed)
    if(dict.found("inletTemperatureList"))
    {
#ifdef OPENFOAMFOUNDATION            
        inletTemperatureList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "inletTemperatureList", dict.subDict("inletTemperatureList")
            )
        );
#elif OPENFOAMESI            
        inletTemperatureList_.reset
        ( 
            Function1<scalar>::New
            (
                "inletTemperatureList", dict.subDict("inletTemperatureList")
            )
        );
#endif        
        
        inletTemperature_ = 
        inletTemperatureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
    else if (dict.found("inletTemperature"))
    {
        inletTemperature_ = readScalar(dict.lookup("inletTemperature"));
    }
    else
    {
        FatalErrorInFunction() 
            << "Coolant channel BC for " << patch().name() << " requires:" << nl
            << "- \"inletTemperatureList\" or" << nl
            << "- \"inletTemperature\"" << abort(FatalError) << endl;
    }

    // Provide pitch over diameter value if Schar model is chosen for Nu
    if ( modelNusselt_ == "Schad" )
    {
        pOverD_ = readScalar(dict.lookup("pitchOverDiameter"));
    }

    if ( modelNusselt_ == "constant" )
    {
        constNusselt_ = readScalar(dict.lookup("nu"));
    }

    if ( modelCP_ == "constant" )
    {
        constCp_ = readScalar(dict.lookup("cp"));
    }

    if ( modelKappa_ == "constant" )
    {
        constKappa_ = readScalar(dict.lookup("k"));
    }

    if ( modelViscosity_ == "constant" )
    {
        constVisc_ = readScalar(dict.lookup("mu"));
    }

    // Assign pointer T coolant if write option enabled by user
    if ( writeCoolantT_ )
    {
        Tcoolant_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "Tcoolant",
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("", dimTemperature, 0),
                "calculated"
            )
        );

        Hcoolant_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "enthalpyCoolant",
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh(),
                dimensionedScalar("", dimEnergy/dimMass, 0),
                "calculated"
            )
        );

        // This patch ID
        label patchID = this->patch().index();

        // Initialize patch field with inletTemperature
        Tcoolant_().boundaryFieldRef()[patchID] = 
            scalarField( this->size(), inletTemperature_);

        // Initialize patch field with inletEnthalpy
        Hcoolant_().boundaryFieldRef()[patchID] = 
            scalarField( this->size(), NaEnthalpy(inletTemperature_));
    }
}

NaCoolantChannelfvPatchScalarField::
NaCoolantChannelfvPatchScalarField
(
    const NaCoolantChannelfvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
#ifdef OPENFOAMFOUNDATION
    fvPatchField<scalar>(ptf, p, iF, mapper, mappingRequired),
#elif OPENFOAMESI
    fvPatchField<scalar>(ptf, p, iF, mapper),
#endif
    kappaName_(ptf.kappaName_),
    T0_(ptf.T0_),
    h_(ptf.h_),
    alpha_(ptf.alpha_),
    pOverD_(ptf.pOverD_),
    A_(ptf.A_),
    relaxHTC_(ptf.relaxHTC_),
    DH_(ptf.DH_),
    massFlowRate_(ptf.massFlowRate_),
    massFlowRateList_(ptf.massFlowRateList_),
    inletTemperature_(ptf.inletTemperature_),
    inletTemperatureList_(ptf.inletTemperatureList_),
    modelKappa_(ptf.modelKappa_),
    constKappa_(ptf.constKappa_),
    modelCP_(ptf.modelCP_),
    constCp_(ptf.constCp_),
    modelViscosity_(ptf.modelViscosity_),
    constVisc_(ptf.constVisc_),
    modelNusselt_(ptf.modelNusselt_),
    constNusselt_(ptf.constNusselt_),
    writeCoolantT_(ptf.writeCoolantT_),
    Tcoolant_(ptf.Tcoolant_),
    heatFlux_(ptf.heatFlux_),
    mapper_(ptf.mapper_)
{}

NaCoolantChannelfvPatchScalarField::
NaCoolantChannelfvPatchScalarField
(
    const NaCoolantChannelfvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(tppsf, iF),
    kappaName_(tppsf.kappaName_),
    T0_(tppsf.T0_),
    h_(tppsf.h_),
    alpha_(tppsf.alpha_),
    pOverD_(tppsf.pOverD_),
    A_(tppsf.A_),
    relaxHTC_(tppsf.relaxHTC_),
    DH_(tppsf.DH_),
    massFlowRate_(tppsf.massFlowRate_),
    massFlowRateList_(),
    inletTemperature_(tppsf.inletTemperature_),
    inletTemperatureList_(),
    modelKappa_(tppsf.modelKappa_),
    constKappa_(tppsf.constKappa_),
    modelCP_(tppsf.modelCP_),
    constCp_(tppsf.constCp_),
    modelViscosity_(tppsf.modelViscosity_),
    constVisc_(tppsf.constVisc_),
    modelNusselt_(tppsf.modelNusselt_),
    constNusselt_(tppsf.constNusselt_),
    writeCoolantT_(tppsf.writeCoolantT_),
    Tcoolant_(),
    heatFlux_(tppsf.heatFlux_),
    mapper_(tppsf.mapper_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void NaCoolantChannelfvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Update mass flow rate
    if (massFlowRateList_.valid())
    {
        massFlowRate_ = 
        massFlowRateList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // Update inlet temperature
    if (inletTemperatureList_.valid())
    {
        inletTemperature_ = 
        inletTemperatureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    if(mapper_ == nullptr)
    {
        mapper_ = &patch().boundaryMesh().mesh().lookupObject<sliceMapper>("sliceMapper");
    }

    if (mapper_->nSlices() <= 0)
    {
        FatalErrorInFunction()
            << "Number of slices cannot be less than 1." << nl
            << "Check that a sliceMapper different than \"none\" is used "
            << "or check the slice definition."
            << abort(FatalError);
    }

    // Compute slice heights
    scalarField sliceHeights = computeSliceHeights();

    // Compute face radii and axial elevation of the face centers
    vectorField centerOfPatchFaces(this->patch().Cf());
    scalarField zetas(this->size(), 0.0);
    scalarField radii(this->size(), 0.0);
    forAll(centerOfPatchFaces, i)
    {
        vector point(centerOfPatchFaces[i]);
        zetas[i] = point.z();
        radii[i] = pow(pow(point.x(),2.0) + pow(point.y(),2.0) , 0.5);
    }

    // Initializing wall temperature
    scalarField Twall = *this;

    // Initializing fuel conductivity and wall heat flux
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fvPatch& patch = this->patch();
    const scalarField& deltaCoeffs = patch.deltaCoeffs();
    scalarField alphaN(this->size(), 0.0);
    scalarField alphaP(this->size(), 0.0);

    if (mesh.foundObject<volTensorField>(kappaName_))
    {
        const fvPatchField<tensor>& k = 
            patch.lookupPatchField<volTensorField, tensor>
            (
                kappaName_
            );
        
        vectorField n = (patch.Sf() / patch.magSf());
        alphaP = ((n & k) & n) * deltaCoeffs;
        heatFlux_ = (patchInternalField() - Twall) * alphaP;
    }
    else
    {
        const fvPatchField<scalar>& k = 
            patch.lookupPatchField<volScalarField, scalar>
            (
                kappaName_
            );

        alphaP = k * deltaCoeffs;
        heatFlux_ = (patchInternalField() - Twall) * alphaP;
    }

    scalarField oldH(this->size(), 0.0);
    scalarField Re(this->size(), 0.0);
    scalarField Pr(this->size(), 0.0);

    // Initialize coolant enthalpy
    scalarField coolantEnthalpy(this->size(), 0.0);

    // Compute film temperature
    scalarField Tfilm = (T0_ + Twall) * 0.5; 

    // Calculate coolant enthalpy for each slice
    forAll(zetas, i)
    {
        if (i == 0)
        {
            coolantEnthalpy[i] = NaEnthalpy(inletTemperature_) + 
                heatFlux_[i] * 2.0 * radii[i] * constant::mathematical::pi * 
                sliceHeights[i] / massFlowRate_;

            T0_[i] =  inletTemperature_ + (coolantEnthalpy[i]-NaEnthalpy(inletTemperature_)) / NaCp( 0.5*(Tfilm[i] + inletTemperature_) );
        }
        else
        {
            coolantEnthalpy[i] = coolantEnthalpy[i - 1] + 
                heatFlux_[i] * 2.0 * radii[i] * constant::mathematical::pi * 
                sliceHeights[i] / massFlowRate_;
            
            T0_[i] =  T0_[i-1] + (coolantEnthalpy[i]-coolantEnthalpy[i-1]) / NaCp( 0.5*(Tfilm[i]+Tfilm[i-1]) );
        }

        // Compute liquid sodium temperature from enthalpy
        // T0_[i] =  inletTemperature_ + (coolantEnthalpy[i]-NaEnthalpy(inletTemperature_)) / NaCp( Tfilm[i] );
        
        // Compute Reynolds
        Re[i] = massFlowRate_ / A_ * DH_ / NaVisc( Tfilm[i] ) ;
        
        // Compute Prandtl
        Pr[i] = NaPrandtl( Tfilm[i] );
    }

    // Compute Peclet number
    scalarField Pe = Re * Pr;

    // Initialize Nusselt number array
    scalarField Nu(zetas.size(), 0.0);

    // Calculating HTC
    forAll(zetas, i)
    {
        if (modelNusselt_ == "constant")
        {
            Nu[i] = constNusselt_;
        }
        // Modified Schad Correlation
        // For triangular sub-channels, applicable for Pe <= 1000 
        else if (modelNusselt_ == "Schad")
        {
            if(Pe[i] >= 1000)
            {
                WarningIn("void NaCoolantChannelfvPatchScalarField::updateCoeff") 
                << nl << " Peclet number of " << Pe[i] << " is outside the " 
                << "applicability range for the Schad model." << endl;
            }

            // Set exponential term to 1 for Pe <= 150
            scalar exp = (Pe[i] <= 150)
            ? 1.0
            : pow(Pe[i]/150, 0.3);

            Nu[i] = 4.496 * (-16.15 + 24.96 * pOverD_ - 8.55 * pow(pOverD_, 2.0) ) * exp;
        }
        // Lyons's law
        // Flow with a surrounding circular tube
        // Applicable for Pr <= 1 and 4000 <= Re <= 3.24e6
        else if (modelNusselt_ == "Lyon")
        {
            if(Re[i] <= 4000 or Re[i] >= 3.24e6)
            {
                WarningIn("void NaCoolantChannelfvPatchScalarField::updateCoeff") 
                << nl << " Reynolds number of " << Re[i] << " is outside the " 
                << "applicability range for the Lyon model." << endl;
            }
            if(Pr[i] >= 1)
            {
                WarningIn("void NaCoolantChannelfvPatchScalarField::updateCoeff") 
                << nl << " Prandtl number of " << Pr[i] << " is outside the " 
                << "applicability range for the Lyon model." << endl;
            }

            Nu[i] = 7.0 + 0.025*pow(max(Pe[i], 200), 0.8);
        }
        // Seban and Shimazaki Correlation
        // For heat transfer from a rod to fluid within a surrounding circular 
        // tube with constant rod wall temperature
        // Applicable for 40 <= Pe <= 1150; 0.0057 <= Pr <= 0.0075 and 5850 <= Re <= 1.78e5
        else if (modelNusselt_ == "SebanShimazaki")
        {
            if(Pe[i] < 40 or Pe[i] > 1150)
            {
                WarningIn("void NaCoolantChannelfvPatchScalarField::updateCoeff") 
                << nl << " Peclet number of " << Pe[i] << " is outside the " 
                << "applicability range for the SebanShimazaki model." << endl;
            }
            if(Re[i] < 5850 or Re[i] > 178000)
            {
                WarningIn("void NaCoolantChannelfvPatchScalarField::updateCoeff") 
                << nl << " Reynolds number of " << Re[i] << " is outside the " 
                << "applicability range for the SebanShimazaki model." << endl;
            }
            if(Pr[i] < 0.0057 or Pr[i] > 0.0075)
            {
                WarningIn("void NaCoolantChannelfvPatchScalarField::updateCoeff") 
                << nl << " Prandtl number of " << Pr[i] << " is outside the " 
                << "applicability range for the SebanShimazaki model." << endl;
            }

            Nu[i] = 5.0 + 0.025 * pow(Pe[i], 0.8);
        }
            else if (modelNusselt_ == "DittusBoelter")
        {
            Nu[i] = 0.023 * pow(Re[i], 0.8) * pow(Pr[i], 0.4);
        }
        else
        {
            FatalErrorInFunction()  
            << "Incorrectly specified Nusselt model name." << nl
            << "Available models are :" << nl  
            << "  - Schad" << nl  
            << "  - Lyon" << nl  
            << "  - SebanShimazaki" << nl << nl  
            << "  - DittusBoelter" << nl << nl  
            << abort(FatalError);
        }

        // Compute HTC
        h_[i] = NaConductivity( Tfilm[i] ) * Nu[i] / DH_;
    }

    // Relaxation of the HTC
    if (sum(oldH) > 0.0)
    {
        h_ = oldH * (1.0 - relaxHTC_) + h_ * relaxHTC_;
    }

    // Assign T coolant field (for postProcessing)
    if ( Tcoolant_.valid() and Hcoolant_.valid() )
    {
        label patchID = this->patch().index();
        Tcoolant_().boundaryFieldRef()[patchID] = T0_;
        Hcoolant_().boundaryFieldRef()[patchID] = coolantEnthalpy;
    }

    // Calculating alphas
    forAll(this->patch().Cf(),faceI)
    {
        alphaN[faceI] = h_[faceI];
    }

    alpha_ = alphaP / (alphaP + alphaN);

    fvPatchField<scalar>::updateCoeffs();
}


void NaCoolantChannelfvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<scalar>::operator=(
            alpha_*patchInternalField() + (1.0 - alpha_) * T0_
    );

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar> >
NaCoolantChannelfvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar> >
NaCoolantChannelfvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1.0 - alpha_) * T0_;
}


tmp<Field<scalar> >
NaCoolantChannelfvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (alpha_ - 1.0) * deltaCoeffs;
}


tmp<Field<scalar> >
NaCoolantChannelfvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (1.0 - alpha_) * deltaCoeffs * T0_;
}


void NaCoolantChannelfvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION

    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    writeEntry(os, "hydraulicDiameter", DH_);
    writeEntry(os, "relaxHTC", relaxHTC_);
    writeEntry(os, "flowArea", A_);
    writeEntry(os, "heatFlux", heatFlux_);
    if (inletTemperatureList_.valid())
    {
        os.writeKeyword("inletTemperatureList") << nl;
        os << tab << token::BEGIN_BLOCK << nl;
        inletTemperatureList_->write(os);
        os << tab << token::END_BLOCK << nl;
    }
    else
    {
        writeEntry(os, "inletTemperature", inletTemperature_);
    }
    if (massFlowRateList_.valid())
    {
        os.writeKeyword("massFlowRateList") << nl;
        os << token::BEGIN_BLOCK << nl;
        massFlowRateList_->write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        writeEntry(os, "massFlowRate", massFlowRate_);
    }
    writeEntry(os, "T0"   , T0_);
    writeEntry(os, "h"    , h_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI

    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    os.writeEntry("hydraulicDiameter", DH_);
    os.writeEntry("relaxHTC", relaxHTC_);
    os.writeEntry("flowArea", A_);
    heatFlux_.writeEntry("heatFlux", os);
    if (inletTemperatureList_.valid())
    {
        os.writeKeyword("inletTemperatureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        inletTemperatureList_->writeData(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        os.writeEntry("inletTemperature", inletTemperature_);
    }
    if (massFlowRateList_.valid())
    {
        os.writeKeyword("massFlowRateList") << nl;
        os << token::BEGIN_BLOCK << nl;
        massFlowRateList_->writeData(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        os.writeEntry("massFlowRate", massFlowRate_);
    }
    T0_.writeEntry("T0", os);
    h_.writeEntry("h", os);
    writeEntry("value", os);
#endif  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, NaCoolantChannelfvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
