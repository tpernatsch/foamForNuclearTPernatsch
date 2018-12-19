/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "volFields.H"
#include "surfaceFields.H"
#include "mappedPatchBase.H"
#include "turbulentFluidThermoModel.H"
#include "mapDistribute.H"
#include "rhoThermo.H" //s

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class solidType>
nusseltThermalBaffle1DFvPatchScalarField<solidType>::
nusseltThermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(p.patch()),
    mixedFvPatchScalarField(p, iF),
    singleDh_(),
    QInfo_(),
    fh_(),
    UName_("U"),
    cellZones1_(),
    cellZones2_(),
    Dh1_(),
    Dh2_(),
    Dh_(),
    lowRe_(),
    highRe_(),
    lNuConst_(),
    lRePrCoeff_(),
    lExpRe_(),
    lExpPr_(),
    tNuConst_(),
    tRePrCoeff_(),
    tExpRe_(),
    tExpPr_(),
    TpReal_(p.size()),
    h_(p.size()),
    Ti_(p.size()),
    fRe_(p.size()),
    Nu_(p.size()),
    TName_("T"),
    baffleActivated_(true),
    thickness_(p.size()),
    solidDict_(),
    solidPtr_()
{}

template<class solidType>
nusseltThermalBaffle1DFvPatchScalarField<solidType>::
nusseltThermalBaffle1DFvPatchScalarField
(
    const nusseltThermalBaffle1DFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedPatchBase(p.patch(), ptf),
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    singleDh_(ptf.singleDh_),
    QInfo_(ptf.QInfo_),
    fh_(ptf.fh_),
    UName_(ptf.UName_),
    cellZones1_(ptf.cellZones1_),
    cellZones2_(ptf.cellZones2_),
    Dh1_(ptf.Dh1_),
    Dh2_(ptf.Dh2_),
    Dh_(ptf.Dh_),
    lowRe_(ptf.lowRe_),
    highRe_(ptf.highRe_),
    lNuConst_(ptf.lNuConst_),
    lRePrCoeff_(ptf.lRePrCoeff_),
    lExpRe_(ptf.lExpRe_),
    lExpPr_(ptf.lExpPr_),
    tNuConst_(ptf.tNuConst_),
    tRePrCoeff_(ptf.tRePrCoeff_),
    tExpRe_(ptf.tExpRe_),
    tExpPr_(ptf.tExpPr_),
    TpReal_(ptf.TpReal_, mapper),
    h_(ptf.h_, mapper),
    Ti_(ptf.Ti_, mapper),
    fRe_(ptf.fRe_, mapper),
    Nu_(ptf.Nu_, mapper),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_, mapper),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_)
{}

template<class solidType>
nusseltThermalBaffle1DFvPatchScalarField<solidType>::
nusseltThermalBaffle1DFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mappedPatchBase(p.patch(), NEARESTPATCHFACE, dict),
    mixedFvPatchScalarField(p, iF),
    singleDh_(dict.lookupOrDefault<bool>("singleDh", false)),
    QInfo_(dict.lookupOrDefault<bool>("QInfo", false)),
    fh_(dict.lookupOrDefault<scalar>("fh", 1.0)),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    cellZones1_(dict.lookupOrDefault<List<word>>("cellZones1", List<word>())),
    cellZones2_(dict.lookupOrDefault<List<word>>("cellZones2", List<word>())),
    Dh1_(dict.lookupOrDefault<scalar>("Dh1", VSMALL)),
    Dh2_(dict.lookupOrDefault<scalar>("Dh2", VSMALL)),
    Dh_(dict.lookupOrDefault<scalar>("Dh", VSMALL)),
    lowRe_(dict.lookupOrDefault<scalar>("lowRe", 1000)),
    highRe_(dict.lookupOrDefault<scalar>("highRe", 2300)),
    lNuConst_(dict.lookupOrDefault<scalar>("lNuConst", 1)),
    lRePrCoeff_(dict.lookupOrDefault<scalar>("lRePrCoeff", 0)),
    lExpRe_(dict.lookupOrDefault<scalar>("lExpRe", 0)),
    lExpPr_(dict.lookupOrDefault<scalar>("lExpPr", 0)),
    tNuConst_(dict.lookupOrDefault<scalar>("tNuConst", 1)),
    tRePrCoeff_(dict.lookupOrDefault<scalar>("tRePrCoeff", 0)),
    tExpRe_(dict.lookupOrDefault<scalar>("tExpRe", 0)),
    tExpPr_(dict.lookupOrDefault<scalar>("tExpPr", 0)),
    TpReal_(),
    h_(),
    Ti_(),
    fRe_(),
    Nu_(),
    TName_("T"),
    baffleActivated_(dict.lookupOrDefault<bool>("baffleActivated", true)),
    thickness_(),
    solidDict_(dict),
    solidPtr_()
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    //- The token issue came from the fact that I was reading TpReal_
    //  by first checking the presence of the "value" field to read
    //  from with dict.found("value"). This somehow messes up the 
    //  subsequent reading of the "value" field for the real purpose
    //  of reading T. Solved it by giving TpReal its own field realValue
    //  in the dict. However, assigning the "value" field to TpReal if no
    //  realValue field is present still works, it is only using
    //  dict.found("value") that messes everything up
    if (dict.found("realValue"))
    {
        TpReal_ = scalarField("realValue", dict, p.size());
    }
    else
    {
        TpReal_ = scalarField("value", dict, p.size());
    }

    if (dict.found("h"))
    {
        h_ = scalarField("h", dict, p.size());
    }
    else
    {
        //- I can't use computeH() here because it would cause
        //- issues when using utilities such as decomposePar
        h_ = scalarField(p.size(), GREAT);
    }

    Ti_ = patchInternalField();

    fRe_= scalarField(p.size(), scalar(0));

    Nu_ = scalarField(p.size(), scalar(0));

    if (dict.found("thickness"))
    {
        thickness_ = scalarField("thickness", dict, p.size());
    }

    if (dict.found("refValue") && baffleActivated_)
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }
}


template<class solidType>
nusseltThermalBaffle1DFvPatchScalarField<solidType>::
nusseltThermalBaffle1DFvPatchScalarField
(
    const nusseltThermalBaffle1DFvPatchScalarField& ptf
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf),
    singleDh_(ptf.singleDh_),
    QInfo_(ptf.QInfo_),
    fh_(ptf.fh_),
    UName_(ptf.UName_),
    cellZones1_(ptf.cellZones1_),
    cellZones2_(ptf.cellZones2_),
    Dh1_(ptf.Dh1_),
    Dh2_(ptf.Dh2_),
    Dh_(ptf.Dh_),
    lowRe_(ptf.lowRe_),
    highRe_(ptf.highRe_),
    lNuConst_(ptf.lNuConst_),
    lRePrCoeff_(ptf.lRePrCoeff_),
    lExpRe_(ptf.lExpRe_),
    lExpPr_(ptf.lExpPr_),
    tNuConst_(ptf.tNuConst_),
    tRePrCoeff_(ptf.tRePrCoeff_),
    tExpRe_(ptf.tExpRe_),
    tExpPr_(ptf.tExpPr_),
    TpReal_(ptf.TpReal_),
    h_(ptf.h_),
    Ti_(ptf.Ti_),
    fRe_(ptf.fRe_),
    Nu_(ptf.Nu_),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_)
{}

template<class solidType>
nusseltThermalBaffle1DFvPatchScalarField<solidType>::
nusseltThermalBaffle1DFvPatchScalarField
(
    const nusseltThermalBaffle1DFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mappedPatchBase(ptf.patch().patch(), ptf),
    mixedFvPatchScalarField(ptf, iF),
    singleDh_(ptf.singleDh_),
    QInfo_(ptf.QInfo_),
    fh_(ptf.fh_),
    UName_(ptf.UName_),
    cellZones1_(ptf.cellZones1_),
    cellZones2_(ptf.cellZones2_),
    Dh1_(ptf.Dh1_),
    Dh2_(ptf.Dh2_),
    Dh_(ptf.Dh_),
    lowRe_(ptf.lowRe_),
    highRe_(ptf.highRe_),
    lNuConst_(ptf.lNuConst_),
    lRePrCoeff_(ptf.lRePrCoeff_),
    lExpRe_(ptf.lExpRe_),
    lExpPr_(ptf.lExpPr_),
    tNuConst_(ptf.tNuConst_),
    tRePrCoeff_(ptf.tRePrCoeff_),
    tExpRe_(ptf.tExpRe_),
    tExpPr_(ptf.tExpPr_),
    TpReal_(ptf.TpReal_),
    h_(ptf.h_),
    Ti_(ptf.Ti_),
    fRe_(ptf.fRe_),
    Nu_(ptf.Nu_),
    TName_(ptf.TName_),
    baffleActivated_(ptf.baffleActivated_),
    thickness_(ptf.thickness_),
    solidDict_(ptf.solidDict_),
    solidPtr_(ptf.solidPtr_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class solidType>
bool nusseltThermalBaffle1DFvPatchScalarField<solidType>::owner() const
{
    const label patchi = patch().index();

    const label nbrPatchi = samplePolyPatch().index();

    return (patchi < nbrPatchi);
}

template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::setProperties()
{
    if (patch().size() != 0)
    {
        const fvMesh& mesh = patch().boundaryMesh().mesh();
        const label cellZoneI = mesh.cellZones().whichZone(patch().faceCells()[0]);

        if (this->owner())
        {   
            if (!singleDh_)
            {
                List<label> cellZones1I(cellZones1_.size());
                forAll(cellZones1_, i)
                {
                    cellZones1I[i] = mesh.cellZones().findZoneID(cellZones1_[i]);
                }
                
                List<label> cellZones2I(cellZones2_.size());
                forAll(cellZones2_, i)
                {
                    cellZones2I[i] = mesh.cellZones().findZoneID(cellZones2_[i]);
                }
                
                forAll(cellZones1I, i)
                {
                    if (cellZoneI == cellZones1I[i])
                    {
                        this->Dh_ = Dh1_;
                        break;
                    }
                }
                
                forAll(cellZones2I, i)
                {
                    if (cellZoneI == cellZones2I[i])
                    {
                        this->Dh_ = Dh2_;
                        break;
                    }
                }

                if (this->Dh_ == VSMALL)
                {
                    FatalErrorIn("scalar nusseltThermalBaffle1DFvPatchScalarField<solidType>::Dh()") <<
                    "Could not find any cell zone with any of the specified names: " << endl <<
                    cellZones1_ << endl << " or : " << cellZones2_ << endl <<exit(FatalError);
                }
            }         
        }
        else
        {
            const fvPatch& nbrPatch =
                patch().boundaryMesh()[samplePolyPatch().index()];

            const nusseltThermalBaffle1DFvPatchScalarField& nbrField =
            refCast<const nusseltThermalBaffle1DFvPatchScalarField>
            (
                nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
            );

            this->singleDh_ = nbrField.singleDh_;

            if (!singleDh_)
            {
                List<label> cellZones1I(nbrField.cellZones1_.size());
                forAll(nbrField.cellZones1_, i)
                {
                    cellZones1I[i] = mesh.cellZones().findZoneID(nbrField.cellZones1_[i]);
                }
                
                List<label> cellZones2I(nbrField.cellZones2_.size());
                forAll(nbrField.cellZones2_, i)
                {
                    cellZones2I[i] = mesh.cellZones().findZoneID(nbrField.cellZones2_[i]);
                }
                
                forAll(cellZones1I, i)
                {
                    if (cellZoneI == cellZones1I[i])
                    {
                        this->Dh_ = nbrField.Dh1_;
                        break;
                    }
                }
                
                forAll(cellZones2I, i)
                {
                    if (cellZoneI == cellZones2I[i])
                    {
                        this->Dh_ = nbrField.Dh2_;
                        break;
                    }
                }

                if (this->Dh_ == VSMALL)
                {
                    FatalErrorIn("scalar nusseltThermalBaffle1DFvPatchScalarField<solidType>::Dh()") <<
                    "Could not find any cell zone with any of the specified names: " << endl <<
                    nbrField.cellZones1_ << endl << " or : " << nbrField.cellZones2_ << endl <<exit(FatalError);
                }
            }
            else
            {
                this->Dh_ = nbrField.Dh_;
            }
            
            this->UName_ = nbrField.UName_;
            this->lowRe_ = nbrField.lowRe_;
            this->highRe_ = nbrField.highRe_;
            this->lNuConst_ = nbrField.lNuConst_;
            this->lRePrCoeff_ = nbrField.lRePrCoeff_;
            this->lExpRe_ = nbrField.lExpRe_;
            this->lExpPr_ = nbrField.lExpPr_;
            this->tNuConst_ = nbrField.tNuConst_;
            this->tRePrCoeff_ = nbrField.tRePrCoeff_;
            this->tExpRe_ = nbrField.tExpRe_;
            this->tExpPr_ = nbrField.tExpPr_;    
        }
    }
}

template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::computeH()
{
    const label patchi = patch().index();
    const volVectorField& U = this->db().objectRegistry::lookupObject<volVectorField>(UName_);
    const rhoThermo& thermo = db().template lookupObject<rhoThermo>("thermophysicalProperties");
    const tmp<volScalarField> nu(thermo.nu()); // Only for the computation of the Reynolds (needs to be the laminar one)
    const tmp<volScalarField> rho(thermo.rho());

    //- Access turbulence model
    const compressible::turbulenceModel& turbModel =
        db().template lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    //- Effective turbulent thermal diffusivity * density
    const tmp<scalarField> alphaEff(turbModel.alphaEff(patchi));

    //- Effective fluid thermal conductivity at the wall
    const scalarField kappaw(turbModel.kappaEff(patchi));

    //- I want some infos on these
    scalarField PrField(this->size(), scalar(0));
    scalarField ReField(this->size(), scalar(0));

    List<label> faceCells(patch().faceCells());
    forAll(faceCells, celli)
    {
        label cell = faceCells[celli];
        scalar Re = mag(U[cell])*Dh_/nu()[cell];
        ReField[celli] = Re;
        scalar Pr = rho()[cell]*nu()[cell]/alphaEff()[celli];
        PrField[celli] = Pr;
        scalar lNu = fh_*lNuConst_ + lRePrCoeff_*pow(Re, lExpRe_)*pow(Pr, lExpPr_);
        scalar tNu = fh_*tNuConst_ + tRePrCoeff_*pow(Re, tExpRe_)*pow(Pr, tExpPr_);
        scalar fRe = fRe_[celli];
        Nu_[celli] = fRe*(tNu)+(1.0-fRe)*(lNu);
        h_[celli] = Nu_[celli]*kappaw[celli]/Dh_;
    }
    /*
    Info << "Patch : " << this->patch().name() << endl
    << " | avg(Dh) : " << Dh_ << " m" << endl
    << " | avg(Re) : " << gAverage(ReField) << endl
    << " | avg(Pr) : " << gAverage(PrField) << endl
    << " | avg(Nu) : " << gAverage(Nu_) << endl
    << " | avg(kE) : " << gAverage(kappaw) << " W/m/K" << endl
    << " | avg(nu) : " << gAverage(nu()) << " m2/s2" << endl
    << " | avg(aE) : " << gAverage(alphaEff()/rho()) << " m2/s2" << endl
    << " | avg(h)  : " << gAverage(h_) << " W/m2/K" << endl;
    */
}

template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::computeFRe()
{
    const volVectorField& U = this->db().objectRegistry::lookupObject<volVectorField>(UName_);
    const rhoThermo& thermo = db().template lookupObject<rhoThermo>("thermophysicalProperties");
    const tmp<volScalarField> nu(thermo.nu()); // (needs to be the laminar one)

    scalar m(1/max((highRe_-lowRe_), SMALL));
    scalar q(-m*lowRe_);

    

    List<label> faceCells(patch().faceCells());
    forAll(faceCells, celli)
    {
        label cell = faceCells[celli];
        scalar Re = mag(U[cell])*Dh_/nu()[cell];
        scalar fRe = m*Re + q;
        fRe_[celli] = max(scalar(0), min(scalar(1), fRe));
    }
}

template<class solidType>
tmp<scalarField> nusseltThermalBaffle1DFvPatchScalarField<solidType>::nbrTpReal() const
{
    const mapDistribute& mapDist = this->mappedPatchBase::map();

    const fvPatch& nbrPatch =
        patch().boundaryMesh()[samplePolyPatch().index()];

    const nusseltThermalBaffle1DFvPatchScalarField& nbrField =
    refCast<const nusseltThermalBaffle1DFvPatchScalarField>
    (
        nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
    );

    tmp<scalarField> tnbrTpReal
    (
        new scalarField(nbrField.TpReal_)
    );

    scalarField& nbrTpReal = tnbrTpReal.ref();
    mapDist.distribute(nbrTpReal);
    return tnbrTpReal;   
}

template<class solidType>
tmp<scalarField> nusseltThermalBaffle1DFvPatchScalarField<solidType>::nbrH() const
{
    const mapDistribute& mapDist = this->mappedPatchBase::map();

    const fvPatch& nbrPatch =
        patch().boundaryMesh()[samplePolyPatch().index()];

    const nusseltThermalBaffle1DFvPatchScalarField& nbrField =
    refCast<const nusseltThermalBaffle1DFvPatchScalarField>
    (
        nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
    );

    tmp<scalarField> tnbrH
    (
        new scalarField(nbrField.h_)
    );

    scalarField& nbrH = tnbrH.ref();
    mapDist.distribute(nbrH);
    return tnbrH; 
}

template<class solidType>
tmp<scalarField> nusseltThermalBaffle1DFvPatchScalarField<solidType>::nbrTi() const
{
    const mapDistribute& mapDist = this->mappedPatchBase::map();

    const fvPatch& nbrPatch =
        patch().boundaryMesh()[samplePolyPatch().index()];

    const nusseltThermalBaffle1DFvPatchScalarField& nbrField =
    refCast<const nusseltThermalBaffle1DFvPatchScalarField>
    (
        nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
    );

    tmp<scalarField> tnbrTi
    (
        new scalarField(nbrField.Ti_)
    );

    scalarField& nbrTi = tnbrTi.ref();
    mapDist.distribute(nbrTi);
    return tnbrTi;   
}

template<class solidType>
const solidType& nusseltThermalBaffle1DFvPatchScalarField<solidType>::solid() const
{
    if (this->owner())
    {
        if (solidPtr_.empty())
        {
            solidPtr_.reset(new solidType(solidDict_));
        }
        return solidPtr_();
    }
    else
    {
        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];

        const nusseltThermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const nusseltThermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        return nbrField.solid();
    }
}

template<class solidType>
tmp<scalarField> nusseltThermalBaffle1DFvPatchScalarField<solidType>::
baffleThickness() const
{
    if (this->owner())
    {
        if (thickness_.size() != patch().size())
        {
            FatalIOErrorInFunction
            (
                solidDict_
            )<< " Field thickness has not been specified "
            << " for patch " << this->patch().name()
            << exit(FatalIOError);
        }

        return thickness_;
    }
    else
    {
        const mapDistribute& mapDist = this->mappedPatchBase::map();

        const fvPatch& nbrPatch =
            patch().boundaryMesh()[samplePolyPatch().index()];
        const nusseltThermalBaffle1DFvPatchScalarField& nbrField =
        refCast<const nusseltThermalBaffle1DFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(TName_)
        );

        tmp<scalarField> tthickness
        (
            new scalarField(nbrField.baffleThickness())
        );
        scalarField& thickness = tthickness.ref();
        mapDist.distribute(thickness);
        return tthickness;
    }
}

template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    if (this->owner())
    {
        thickness_.autoMap(m);
        TpReal_.autoMap(m);
        h_.autoMap(m);
        Ti_.autoMap(m);
    }
}


template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const nusseltThermalBaffle1DFvPatchScalarField& tiptf =
        refCast<const nusseltThermalBaffle1DFvPatchScalarField>(ptf);

    if (this->owner())
    {
        thickness_.rmap(tiptf.thickness_, addr);
        TpReal_.rmap(tiptf.TpReal_, addr);
        h_.rmap(tiptf.h_, addr);
        Ti_.rmap(tiptf.Ti_, addr);
    }
}

template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we are inside initEvaluate/evaluate there might be 
    // processor comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    if (baffleActivated_)
    {
        //- Set baffle properties (hydraulic diameter, correlation parameters)
        //  This is done here, only at the first time step, rather than in the 
        //  constructor because in the constructor it causes issues when
        //  using decomposePar
        if (this->db().time().timeIndex() - this->db().time().startTimeIndex() == 1)
        {
            setProperties();
        }

        //- Compute weighting factors between laminar and turbulent correlations
        computeFRe();

        //- Update the heat transfer coefficient (the one with the Nusselt correlation)
        computeH();

        //- Update internal temperature field (class member)
        Ti_ = patchInternalField();

        //- Access turbulence model and compute fluid heat transfer
        //  coefficient at the wall (as seen by the solver)
        const compressible::turbulenceModel& turbModel =
            db().template lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );  
        const scalarField kappaw(turbModel.kappaEff(patch().index()));
        const scalarField kDelta(patch().deltaCoeffs()*kappaw);

        //- Compute baffle heat transfer coefficient (KDeltaSolid)
        scalarField kappas(patch().size(), 0.0);
        tmp<scalarField> avgTBaffle((TpReal_ + nbrTpReal())/2);
        forAll(kappas, i)
        {
            kappas[i] = solid().kappa(0.0, (avgTBaffle()[i]));
        }
        scalarField KDeltaSolid(kappas/baffleThickness());        
    
        //- Combined heat transfer coefficient for baffle + (fluid on neighbour side)
        scalarField H
        (
            1/
            (
                (1/max(KDeltaSolid, SMALL))
                +
                (1/max(nbrH(), SMALL))
            )
        );

        //- Compute real patch temperature
        tmp<scalarField> vF(H/(H + h_));
        TpReal_ = vF()*nbrTi()+(1-vF())*Ti_;

        //- The following variables will define the actual patch temperature 
        //  value as : 
        //  Tp = refValue()*valueFraction() + 
        //     + (1-valueFraction())*(this->patchInternalField() + 
        //     + refGrad()*patch().deltaCoeffs())
        valueFraction() = 0;
        refValue() = 0;
        refGrad() = (h_/kappaw)*(TpReal_-Ti_);

        if (QInfo_)
        {
            Info    << "Baffle " << patch().name()
                    << ", avg(Q) = " << gAverage(KDeltaSolid*(TpReal_-nbrTpReal())) << " W/m2" << endl;
        }

        //- More infos if debugging (cout rather than Info to 
        //  see individual core info)
        if (debug)
        {
            scalarField updatedTp((h_/kDelta)*(TpReal_-Ti_)+Ti_);
            scalar Qf = gAverage(kDelta*(Ti_-updatedTp));
            scalar Qfh = gAverage(h_*(Ti_-TpReal_));
            scalar Qs = gAverage(KDeltaSolid*(TpReal_-nbrTpReal()));
            scalar nbrQ = gAverage(H*(TpReal_-nbrTi()));
            scalar nbrQfh = gAverage(nbrH()*(nbrTpReal()-nbrTi()));
            cout << "[" << Pstream::myProcNo() << "] | Patch name and size          : " << patch().name() << " - " << patch().size() << "\n";
            cout << "[" << Pstream::myProcNo() << "] | Dh : " << this->Dh_ << " m " << endl; //- NuConst : " << NuConst_ << " - RePrCoeff : " << RePrCoeff_ << " - expRe : " << expRe_ << " - expPr : " << expPr_ << "\n";
            cout << "[" << Pstream::myProcNo() << "] | Avg patch T                  : " << gAverage(updatedTp) << " K\n" << endl;
            cout << "[" << Pstream::myProcNo() << "] | Avg patch real T             : " << gAverage(TpReal_) << " K\n" << endl;
            cout << "[" << Pstream::myProcNo() << "] | Avg Q fluid (via kDelta)     : " << Qf << " W/m2\n" << endl;
            cout << "[" << Pstream::myProcNo() << "] | Avg Q fluid (via h(Nu))      : " << Qfh << " W/m2\n" << endl;
            cout << "[" << Pstream::myProcNo() << "] | Avg Q baffle                 : " << Qs << " W/m2\n" << endl;
            cout << "[" << Pstream::myProcNo() << "] | Avg nbrQ (baffle + nbrH(Nu)) : " << nbrQ << " W/m2\n" << endl;
            cout << "[" << Pstream::myProcNo() << "] | Avg nbrQ fluid (via h(Nu))   : " << nbrQfh << " W/m2\n" << endl;
            cout << "\n";
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    mixedFvPatchScalarField::updateCoeffs();
}

template<class solidType>
void nusseltThermalBaffle1DFvPatchScalarField<solidType>::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);

    TpReal_.writeEntry("realValue", os);

    os.writeKeyword("QInfo")<< QInfo_
            << token::END_STATEMENT << nl;

    os.writeKeyword("fh")<< fh_
            << token::END_STATEMENT << nl;

    h_.writeEntry("h", os);

    mappedPatchBase::write(os);

    if (this->owner())
    {
        baffleThickness()().writeEntry("thickness", os);
        solid().write(os);
        os.writeKeyword("singleDh")<< singleDh_ 
        << token::END_STATEMENT << nl;
        if (!singleDh_)
        {
            os.writeKeyword("cellZones1")<< cellZones1_
            << token::END_STATEMENT << nl;
            os.writeKeyword("Dh1")<< Dh1_
                << token::END_STATEMENT << nl;
            os.writeKeyword("cellZones2")<< cellZones2_
                << token::END_STATEMENT << nl;
            os.writeKeyword("Dh2")<< Dh2_
                << token::END_STATEMENT << nl;
        }
    }

    if (Dh_ != VSMALL)
    {
        os.writeKeyword("Dh")<< Dh_
        << token::END_STATEMENT << nl;
    } 
    
    if (this->owner())
    {
        os.writeKeyword("lowRe")<< lowRe_
            << token::END_STATEMENT << nl;
        os.writeKeyword("highRe")<< highRe_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lNuConst")<< lNuConst_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lRePrCoeff")<< lRePrCoeff_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lExpRe")<< lExpRe_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lExpPr")<< lExpPr_
            << token::END_STATEMENT << nl;
        os.writeKeyword("tNuConst")<< tNuConst_
            << token::END_STATEMENT << nl;
        os.writeKeyword("tRePrCoeff")<< tRePrCoeff_
            << token::END_STATEMENT << nl;
        os.writeKeyword("tExpRe")<< tExpRe_
            << token::END_STATEMENT << nl;
        os.writeKeyword("tExpPr")<< tExpPr_
            << token::END_STATEMENT << nl;
        os.writeKeyword("UName")<< UName_
            << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
