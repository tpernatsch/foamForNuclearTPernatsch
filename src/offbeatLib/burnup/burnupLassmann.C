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

    OpenFOAM is distributed in the hope that it will be useful, BUt WITHOUT
    ANY WARRANTY+ without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "burnupLassmann.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "Tuple2.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(burnupLassmann, 0);
    addToRunTimeSelectionTable
    (
        burnup, 
        burnupLassmann, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::burnupLassmann::setXSForIntEnrichment
(
    const int eRounded,
    scalarField& sf,
    scalarField& sc
) const
{
    // assume sf, sc already sized to Nelements_ and initialised to 0

    switch (eRounded)
    {
        case 2:
        // enr : U235/U: 2 wt-%
        // Fission
        sf[iU235_]  = 50.3954e-24;
        sf[iU236_]  = 0.264092e-24;
        sf[iU238_]  = 0.0867702e-24;
        sf[iNp237_] = 0.494232e-24;
        sf[iPu238_] = 2.52238e-24;
        sf[iPu239_] = 119.237e-24;
        sf[iPu240_] = 0.533673e-24;
        sf[iPu241_] = 127.444e-24;
        sf[iPu242_] = 0.393512e-24;
        sf[iAm241_] = 1.10784e-24;
        sf[iAm243_] = 0.398466e-24;
        sf[iCm242_] = 0.538437e-24;
        sf[iCm243_] = 88.8823e-24;
        sf[iCm244_] = 0.798091e-24;

        // Capture :
        sc[iU235_]  = 11.1121e-24;
        sc[iU236_]  = 6.49929e-24;
        sc[iU238_]  = 0.916154e-24;
        sc[iNp237_] = 34.9828e-24;
        sc[iPu238_] = 39.3757e-24;
        sc[iPu239_] = 65.6719e-24;
        sc[iPu240_] = 115.326e-24; // <--- TODO we miss bu evolution wrt TUBRNP
        sc[iPu241_] = 45.8265e-24;
        sc[iPu242_] = 25.0506e-24;
        sc[iAm241_] = 115.076e-24;
        sc[iAm243_] = 48.3073e-24;
        sc[iCm242_] = 4.45425e-24;
        sc[iCm243_] = 14.9071e-24;
        sc[iCm244_] = 17.8943e-24;

        break;

        case 3:
        // enr : U235/U: 3 wt-%
        // Fission :
        sf[iU235_]  = 43.7224e-24;
        sf[iU236_]  = 0.264443e-24;
        sf[iU238_]  = 0.0890755e-24;
        sf[iNp237_] = 0.507123e-24;
        sf[iPu238_] = 2.38196e-24;
        sf[iPu239_] = 104.189e-24;
        sf[iPu240_] = 0.546862e-24;
        sf[iPu241_] = 110.574e-24;
        sf[iPu242_] = 0.404023e-24;
        sf[iAm241_] = 1.06063e-24;
        sf[iAm243_] = 0.407574e-24;
        sf[iCm242_] = 0.510015e-24;
        sf[iCm243_] = 81.4531e-24;
        sf[iCm244_] = 0.80456e-24;

        // Capture :
        sc[iU235_]  = 9.97339e-24;
        sc[iU236_]  = 6.11708e-24;
        sc[iU238_]  = 0.896241e-24;
        sc[iNp237_] = 33.3114e-24;
        sc[iPu238_] = 33.6815e-24;
        sc[iPu239_] = 57.7786e-24;
        sc[iPu240_] = 113.395e-24; // <--- TODO we miss bu evolution wrt TUBRNP
        sc[iPu241_] = 39.8099e-24;
        sc[iPu242_] = 26.8681e-24;
        sc[iAm241_] = 104.742e-24;
        sc[iAm243_] = 48.4411e-24;
        sc[iCm242_] = 4.30791e-24;
        sc[iCm243_] = 13.448e-24;
        sc[iCm244_] = 18.1368e-24;

        break;

        case 4:
        // enr : U235/U: 4 wt-%
        // Fission :

        sf[iU235_]  = 37.4386e-24;
        sf[iU236_]  = 0.265747e-24;
        sf[iU238_]  = 0.0912251e-24;
        sf[iNp237_] = 0.51919e-24;        
        sf[iPu238_] = 2.346e-24;
        sf[iPu239_] = 89.736e-24;
        sf[iPu240_] = 0.611e-24;
        sf[iPu241_] = 96.887e-24;
        sf[iPu242_] = 0.475e-24;
        sf[iAm241_] = 1.01663e-24;
        sf[iAm243_] = 0.415918e-24;
        sf[iCm242_] = 0.483447e-24;
        sf[iCm243_] = 74.4952e-24;
        sf[iCm244_] = 0.80861e-24;

        // Capture :
        sc[iU235_]  = 8.8814e-24;
        sc[iU236_]  = 5.84883e-24;
        sc[iU238_]  = 0.876741e-24;
        sc[iNp237_] = 31.5026e-24;
        sc[iPu238_] = 29.042e-24;
        sc[iPu239_] = 49.923e-24;
        sc[iPu240_] = 110.857e-24;
        sc[iPu241_] = 31.682e-24;
        sc[iPu242_] = 33.693e-24;
        sc[iAm241_] = 95.0178e-24;
        sc[iAm243_] = 48.349e-24;
        sc[iCm242_] = 4.19775e-24;
        sc[iCm243_] = 12.0864e-24;
        sc[iCm244_] = 18.3196e-24;

        break;
        
        case 5:

        // enr : U235/U: 5 wt-%
        // Fission :
        sf[iU235_]  = 32.4072e-24;
        sf[iU236_]  = 0.266853e-24;
        sf[iU238_]  = 0.0934406e-24;
        sf[iNp237_] = 0.531233e-24;
        sf[iPu238_] = 2.251e-24;
        sf[iPu239_] = 80.068e-24;
        sf[iPu240_] = 0.622e-24;
        sf[iPu241_] = 85.617e-24;
        sf[iPu242_] = 0.483e-24;
        sf[iAm241_] = 0.980386e-24;
        sf[iAm243_] = 0.424378e-24;
        sf[iCm242_] = 0.464009e-24;
        sf[iCm243_] = 68.6188e-24;
        sf[iCm244_] = 0.814002e-24;

        // Capture :
        sc[iU235_]  = 7.98174e-24;
        sc[iU236_]  = 5.56608e-24;
        sc[iU238_]  = 0.860247e-24;
        sc[iNp237_] = 29.8341e-24;
        sc[iPu238_] = 25.039e-24;
        sc[iPu239_] = 45.034e-24;
        sc[iPu240_] = 108.836e-24;
        sc[iPu241_] = 27.836e-24;
        sc[iPu242_] = 33.799e-24;
        sc[iAm241_] = 86.6297e-24;
        sc[iAm243_] = 47.9882e-24;
        sc[iCm242_] = 4.06005e-24;
        sc[iCm243_] = 10.9582e-24;
        sc[iCm244_] = 18.0668e-24;

        break;

        default:
            FatalErrorInFunction
                << "Unsupported enrichment (rounded) " << eRounded
                << " – expected 2, 3, 4 or 5" << abort(FatalError);
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burnupLassmann::burnupLassmann
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& burnupDict
)
:
    burnupFromPower(mesh, mat, burnupDict),
    formFactor_
    (  
        IOobject
        (
            "formFactor", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::READ_IF_PRESENT, 
            IOobject::AUTO_WRITE
        ), 
        mesh_, 
        dimensionedScalar("formFactor"   , dimless, 1),    
        zeroGradientFvPatchField<scalar>::typeName
    ),
    BukgU_
    ( 
        IOobject
        (
            "Bu_kgU", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::READ_IF_PRESENT, 
            IOobject::AUTO_WRITE
        ), 
        mesh_, 
        dimensionedScalar("Bu_kgU"   , dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    Nelements_(14),

    // ******************* FIELDS TO TRACK COMPOSITIONS *******************   //   

    // Local atomic density
    N_U235_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_U235", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),  
    N_U236_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_U236", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),  
    N_Np237_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Np237", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),   
    N_Pu238_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Pu238", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),  
    N_U238_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_U238", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Pu239_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Pu239", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Pu240_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Pu240", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Pu241_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Pu241", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Pu242_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Pu242", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Am241_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Am241", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Am243_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Am243", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Cm242_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Cm242", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Cm243_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Cm243", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    N_Cm244_
    (
        createOrLookup<scalar>
        (
            mesh_, 
            "N_Cm244", 
            dimless/dimVolume, 
            0.0,
            "zeroGradient"
        )
    ),
    iU235_ (0),
    iU236_ (1),
    iNp237_(2),
    iPu238_(3),
    iU238_ (4),
    iPu239_(5),
    iPu240_(6),
    iPu241_(7),
    iPu242_(8),
    iAm241_(9),
    iAm243_(10),
    iCm242_(11),
    iCm243_(12),
    iCm244_(13),
    sf_(Nelements_,0.0),
    sc_(Nelements_,0.0),
    sa_th_(Nelements_,0.0),
    sn2n_(0.0),
    rIn_(mesh_.nCells(),GREAT),
    rOut_(mesh_.nCells(),-GREAT),
    enrichment_(mesh_.nCells(), 0.0),
    densityFraction_(mesh_.nCells(), 0.0),
    mapper_(nullptr)
{   

    precision_ = burnupDict_.lookupOrDefault("convergencePrecision", 1e-2);

    const labelListList& matAddrList(mat_.matAddrList());

    // Fuel materials are supposed to have information about enrichement and 
    // density fractions
    forAll(mat_.materialsList(), i)
    {   
        const materialModel& materialI(mat_.materialsList()[i]);
        const word materialName(materialI.name());
        const labelList& addr(matAddrList[i]);
        const dictionary& materialDict(materialI.materialModelDict());

        bool isFuel(isA<fuelMaterial>(mat_.materialsList()[i]));

        if(isFuel)
        {
            const UO2& fuelMat = 
                refCast<const UO2>(mat_.materialsList()[i]);

            scalar enrichment(fuelMat.enrichment());
            scalar densityFrac(fuelMat.densityFrac());

            // At the first step : calculate compositions
            
            const scalar weightU = enrichment*235.044 + (1-enrichment)*238.051;
            const scalar weightUO2 = 2*15.999 + weightU;
            const scalar N_tot_ini = 10960e-3*densityFrac/weightUO2*6.022e23; 
            const scalar N5 = N_tot_ini*enrichment;
            const scalar N8 = N_tot_ini-N5;

            // Initialize local concentration (only if not set before, i.e. only 
            // if their current cell value is > 0)
            forAll(addr, i)
            {
                const label cellI = addr[i];
                
                if(N_U235_[cellI] < SMALL)
                {
                    N_U235_[cellI] = N5; 
                }
                
                if(N_U238_[cellI] < SMALL)
                {
                    N_U238_[cellI] = N8;
                }

                enrichment_[cellI] = enrichment;
                densityFraction_[cellI] = densityFrac;
            }  
        }
             
        // Compute inner and outer radii
        calcRinRout(addr);

        // Initialize nuclear data
        initNuclearData(addr);
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burnupLassmann::~burnupLassmann()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::burnupLassmann::correct()
{
    if(mapper_ == nullptr)
    {
        mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
    }

    // Update average burnup
    burnupFromPower::correct();

    // Update inner/outer slice dimensions in case of moving mesh
    if(mesh_.foundObject<fvMesh>("referenceMesh"))
    {
        const labelListList& matAddrList(mat_.matAddrList());
        forAll(mat_.materialsList(), i)
        {   
            const labelList& addr(matAddrList[i]);

            // Compute inner and outer radii
            calcRinRout(addr);
        }
    }

    // Number of cells
    scalar Ncells = mesh_.nCells();

    // Get number of slices
    const scalar nSlices(mapper_->nSlices());

    scalarField Bu_av(nSlices, 0.0);
    scalarField Bu_av_old(nSlices, 0.0);
    scalarField volPerSlice(nSlices, 0.0);

    // Burnup in MWd/MTUO2
    scalarField& Bui = Bu_.ref();
    const scalarField& BuOldi = Bu_.oldTime().internalField();

    // Reference to volumes
    const scalarField& V(mesh_.V());

    const PtrList<labelList>& sliceAddr(mapper_->sliceAddrList());

    forAll(sliceAddr, sliceID)
    {  
        const labelList& addr(sliceAddr[sliceID]);

        scalarField Vslice(V, addr);
        scalarField Buslice(Bui, addr);
        scalarField BuOldislice(BuOldi, addr);

        Bu_av[sliceID] = gSum(Buslice*Vslice)/gSum(Vslice)/0.8815;
        Bu_av_old[sliceID] = gSum(BuOldislice*Vslice)/gSum(Vslice)/0.8815;
        volPerSlice[sliceID] = gSum(Vslice);
    }

    // Prepare burnup average increment field
    scalarField dBu(Ncells, 0.0);
    forAll(sliceAddr, sliceID)
    {
        const labelList& addr(sliceAddr[sliceID]);

        forAll(addr, addrI)
        {
            const label cellI(addr[addrI]);
           
            dBu[cellI] =
            (Bu_av[sliceID] - Bu_av_old[sliceID])/1e3;
        }
    }

    scalar averageBu(0);

// READ INPUT AND OLD TIME CONCENTRATIONS 

    const List<bool>& isFuelList(mapper_->isFuel());
    forAll(sliceAddr, sliceID)
    {
        if(isFuelList[sliceID])
        {
            //TOD0 enrichment and density fraction initials or current?
            scalar enrichmentAverage(0.0);
            scalar densityFractionAverage(0.0);

            const List<label>& addrSlice_(sliceAddr[sliceID]);
            scalar localSize(addrSlice_.size());

            scalarField V_local(V, addrSlice_);

            averageBu *= 0;
            forAll(addrSlice_, i)
            {
                const label cellI = addrSlice_[i];
                averageBu += Bui[cellI]*V[cellI]; // cm³  
                enrichmentAverage += enrichment_[cellI]*V[cellI];
                densityFractionAverage += densityFraction_[cellI]*V[cellI];
            }

            reduce(averageBu, sumOp<scalar>());
            reduce(enrichmentAverage, sumOp<scalar>());
            reduce(densityFractionAverage, sumOp<scalar>());

            averageBu /= volPerSlice[sliceID];
            enrichmentAverage /= volPerSlice[sliceID];
            densityFractionAverage /= volPerSlice[sliceID];

            // Old concentrations and total volume
            scalarRectangularMatrix N_Old_(Nelements_,localSize);
            scalarField N_average_Old_(Nelements_,0.0);
            scalarField N_average_(Nelements_,0.0);

            forAll(addrSlice_, i)
            {
                const label cellI = addrSlice_[i];

                N_Old_[iU235_][i]  = N_U235_.oldTime().internalField()[cellI];
                N_Old_[iU236_][i]  = N_U236_.oldTime().internalField()[cellI];
                N_Old_[iNp237_][i] = N_Np237_.oldTime().internalField()[cellI];
                N_Old_[iPu238_][i] = N_Pu238_.oldTime().internalField()[cellI];
                N_Old_[iU238_][i]  = N_U238_.oldTime().internalField()[cellI];
                N_Old_[iPu239_][i] = N_Pu239_.oldTime().internalField()[cellI];
                N_Old_[iPu240_][i] = N_Pu240_.oldTime().internalField()[cellI];
                N_Old_[iPu241_][i] = N_Pu241_.oldTime().internalField()[cellI];
                N_Old_[iPu242_][i] = N_Pu242_.oldTime().internalField()[cellI];
                N_Old_[iAm241_][i] = N_Am241_.oldTime().internalField()[cellI];
                N_Old_[iAm243_][i] = N_Am243_.oldTime().internalField()[cellI];
                N_Old_[iCm242_][i] = N_Cm242_.oldTime().internalField()[cellI];
                N_Old_[iCm243_][i] = N_Cm243_.oldTime().internalField()[cellI];
                N_Old_[iCm244_][i] = N_Cm244_.oldTime().internalField()[cellI];
                
                for (int isotope = 0; isotope < Nelements_; isotope ++)
                {
                    N_average_Old_[isotope] += N_Old_[isotope][i]*V_local[i];
                }

            }

            reduce(N_average_Old_, sumOp<scalarField>());

            N_average_Old_ /= volPerSlice[sliceID];

            N_average_ = N_average_Old_;

    // NUCLEAR DATA

            /*
                TODO : 
                Here the xs of Pu-240 (Bu-dependent) should be corrected !
            */

            // Compute absorption xs
            scalarField sa_(sc_+sf_);

            // Decay constants
            scalarField lambda_(Nelements_,0.0);
            lambda_[iPu238_] = log(2.) / (8.77e+01  * 365.);
            lambda_[iPu241_] = log(2.) / (1.43e+01  * 365.);
            lambda_[iPu242_] = log(2.) / (3.73e+05  * 365.);
            lambda_[iAm241_] = log(2.) / (4.32e+02  * 365.);
            lambda_[iAm243_] = log(2.) / (7.37e+03  * 365.);
            lambda_[iCm242_] = log(2.) / (1.63e+02        );
            lambda_[iCm244_] = log(2.) / (1.81e+01  * 365.);

            // Branching ratios for the Am242
            const scalar BR_EC_Am242 = 0.173;
            const scalar BR_ng_Am242 = 0.88;
            const scalar BR_beta_Am242 = 0.827;
            const scalar BRPu2 = BR_EC_Am242   * BR_ng_Am242; 
            const scalar BRCm2 = BR_beta_Am242 * BR_ng_Am242;

    // RADIAL SHAPE FUNCTIONS

            scalarField frU_(localSize, 0.0);
            scalarField frPu_(localSize, 0.0);  

            scalar normU_ = 0;
            scalar normPu_ = 0;

            const volVectorField& C(mesh_.C());

            // Calculate volume-averaged radial shape functions
            forAll(addrSlice_, i)
            {
                const label cellI = addrSlice_[i];

                scalar r = pow(pow(C[cellI][0],2)+pow(C[cellI][1],2),0.5);

                frU_[i]  = (1 + 3.45*exp(-3.0*pow((rOut_[cellI] - r)*1e3, 0.45 )));
                frPu_[i] = (1 + 0.94*exp(-3.0*pow((rOut_[cellI] - r)*1e3, 0.45 )));

                // Calculate average frU and frPu
                normU_ += frU_[i]*V_local[i];
                normPu_ += frPu_[i]*V_local[i];
            }

            reduce(normU_, sumOp<scalar>());
            reduce(normPu_, sumOp<scalar>());

            normU_ /= volPerSlice[sliceID];
            normPu_ /= volPerSlice[sliceID];

            frU_ = frU_/normU_;
            frPu_ = frPu_/normPu_;

    // CONVERGENCE LOOP      
            
            scalarRectangularMatrix N_old_iteration_(Nelements_, localSize);
            scalarRectangularMatrix N_(Nelements_, localSize);
            scalarField Qprofile(formFactor_, addrSlice_);  
            N_ = N_Old_;

            scalarField dBu_local_no_profile(dBu, addrSlice_);
            scalarField dBu_local_profile = dBu_local_no_profile*Qprofile;

            const scalar alpha = 3.8e-22; // conversion factor

            scalar E_fis = 3.6156e-22; // MWd
            scalar SigmaF_tot_ = 0;
            scalar SigmaF_tot_old = 0;
            scalarField SigmaF_(localSize,0.0);

            for (int isotope = 0; isotope < Nelements_; isotope ++)
            {
                SigmaF_tot_ += sf_[isotope] * N_average_[isotope]; // cm^-1
                SigmaF_tot_old += sf_[isotope] * N_average_Old_[isotope]; // cm^-1

                forAll(addrSlice_, i)
                {
                    SigmaF_[i] += sf_[isotope]*N_[isotope][i];
                }
            }

            SigmaF_tot_=SigmaF_tot_old;

            scalar deltaT = mesh_.time().deltaT().value()/24/3600;
            scalar difference = VGREAT;

            while (difference > precision_)
            {
                scalar rho_hm = 10960*densityFractionAverage*0.8815;

                N_old_iteration_ = N_;
                scalarField S(Nelements_,0.0);

                // Matrix creation and new compositions calculation
                forAll(addrSlice_, i)
                {
                    // const label cellI = addrSlice_[i];

                    scalar A = 10960*densityFractionAverage*0.8815/1e6/(alpha*SigmaF_[i]);


                    // Matrix creation
                    scalarSquareMatrix M_(Nelements_,0.0);
                    scalarSquareMatrix Absorption_(Nelements_,0.0);
                    scalarSquareMatrix Capture_(Nelements_,0.0);
                    scalarSquareMatrix Decay_loss_(Nelements_,0.0);
                    scalarSquareMatrix Decay_creation_(Nelements_,0.0);
                    List<double> B_(Nelements_);
                    List<double> X_(Nelements_);

                    for (int isotope = 0; isotope < Nelements_; isotope ++)
                    {
                        Decay_loss_[isotope][isotope]=lambda_[isotope]*deltaT;
                        B_[isotope]=N_Old_[isotope][i];
                        if (isotope == iU238_)
                        {
                            Absorption_[isotope][isotope]=1+(frU_[i]*sa_[isotope]+sn2n_)*A*dBu_local_profile[i];
                        }
                        else if (isotope == iPu240_)
                        {
                            Absorption_[isotope][isotope]=1+frPu_[i]*sa_[isotope]*A*dBu_local_profile[i];
                        }
                        else
                        {
                            Absorption_[isotope][isotope]=1+1*sa_[isotope]*A*dBu_local_profile[i];
                        }
                    }

                    Capture_[ iU236_  ][ iU235_  ] = -1            * sc_[iU235_]  * A * dBu_local_profile[i];
                    Capture_[ iNp237_ ][ iU236_  ] = -1            * sc_[iU236_]  * A * dBu_local_profile[i];
                    Capture_[ iPu238_ ][ iNp237_ ] = -1            * sc_[iNp237_] * A * dBu_local_profile[i];
                    Capture_[ iPu239_ ][ iPu238_ ] = -1            * sc_[iPu238_] * A * dBu_local_profile[i];
                    Capture_[ iPu239_ ][ iU238_  ] = -frU_[i]  * sc_[iU238_]  * A * dBu_local_profile[i];
                    Capture_[ iPu240_ ][ iPu239_ ] = -1            * sc_[iPu239_] * A * dBu_local_profile[i];
                    Capture_[ iPu241_ ][ iPu240_ ] = -frPu_[i] * sc_[iPu240_] * A * dBu_local_profile[i];
                    Capture_[ iPu242_ ][ iPu241_ ] = -1            * sc_[iPu241_] * A * dBu_local_profile[i];
                    Capture_[ iPu242_ ][ iAm241_ ] = -BRPu2        * sc_[iAm241_] * A * dBu_local_profile[i];
                    Capture_[ iAm243_ ][ iPu242_ ] = -1            * sc_[iPu242_] * A * dBu_local_profile[i];
                    Capture_[ iCm242_ ][ iAm241_ ] = -BRCm2        * sc_[iAm241_] * A * dBu_local_profile[i];
                    Capture_[ iCm243_ ][ iCm242_ ] = -1            * sc_[iCm242_] * A * dBu_local_profile[i];
                    Capture_[ iCm244_ ][ iAm243_ ] = -1            * sc_[iAm243_] * A * dBu_local_profile[i];
                    Capture_[ iCm244_ ][ iCm243_ ] = -1            * sc_[iCm243_] * A * dBu_local_profile[i];
                  
                    Capture_[ iNp237_][ iU238_ ] = -1              * sn2n_        * A * dBu_local_profile[i]; // is not really a capture

                    Decay_creation_[iNp237_][iAm241_] = -lambda_[iAm241_] * deltaT;
                    Decay_creation_[iPu238_][iCm242_] = -lambda_[iCm242_] * deltaT;
                    Decay_creation_[iPu240_][iCm244_] = -lambda_[iCm244_] * deltaT;
                    Decay_creation_[iAm241_][iPu241_] = -lambda_[iPu241_] * deltaT;
                    
                    M_=Absorption_ + Capture_ + Decay_loss_ + Decay_creation_;
                    
                    solve(X_,M_,B_);

                    SigmaF_[i]=0;
                    for (int isotope = 0; isotope < Nelements_; isotope ++)
                    {
                        N_[isotope][i] = X_[isotope];
                        S[isotope] += N_[isotope][i]*V_local[i];
                        SigmaF_[i] += sf_[isotope]*N_[isotope][i];
                    }
                }

                reduce(S, sumOp<scalarField>());
                S /= volPerSlice[sliceID];

                N_average_ = max(S, SMALL);
                
                // Flux, burnup and power profiles calculation
                SigmaF_tot_ = 0;
                for (int isotope = 0; isotope < Nelements_; isotope ++)
                {
                    SigmaF_tot_ += sf_[isotope] * N_average_[isotope];
                }

                // Get flux normalized distribution from registry
                scalarField fluxProfile; 
                if( mesh_.foundObject<volScalarField>("neutronFlux0") )
                {
                    fluxProfile = 
                    scalarField(mesh_.lookupObject<volScalarField>("neutronFlux0")(), addrSlice_);
                }
                else
                {
                    FatalErrorIn("Foam::burnupLassmann::correct()")
                    << "This burnup model needs a neutronics subsolver to "
                    << "be set. " << exit(FatalError);
                }

                // Compute power profile and normalize it by its average
                Qprofile = fluxProfile * SigmaF_;
                Qprofile /= gSum(Qprofile * V_local)/volPerSlice[sliceID];

                // Apply Q profile to local bu
                dBu_local_profile = Qprofile*dBu_local_no_profile;

                // Difference calculation for convergence  
                scalar dif(0);
                scalar maxdif(0);
                for (int isotope = 0; isotope < Nelements_; isotope ++)
                {
                    forAll(addrSlice_, i)
                    {
                        // const label cellI = addrSlice_[i];
                        
                        if (N_old_iteration_[isotope][i]>SMALL)
                        {
                            dif=mag((N_[isotope][i]-N_old_iteration_[isotope][i])/N_old_iteration_[isotope][i]);

                            if (dif>maxdif)
                            {
                                maxdif=dif;
                            }
                        }
                    }
                }

                reduce(maxdif, maxOp<scalar>());

                difference = maxdif;
            }

    // COMPOSITIONS WRITING

            forAll(addrSlice_, i)
            {
                const label cellI = addrSlice_[i]; 
                N_U235_[cellI]  =   N_[iU235_][i] ;
                N_U236_[cellI]  =   N_[iU236_][i] ;
                N_Np237_[cellI] =   N_[iNp237_][i];
                N_Pu238_[cellI] =   N_[iPu238_][i];
                N_U238_[cellI]  =   N_[iU238_][i] ;
                N_Pu239_[cellI] =   N_[iPu239_][i];
                N_Pu240_[cellI] =   N_[iPu240_][i];
                N_Pu241_[cellI] =   N_[iPu241_][i];
                N_Pu242_[cellI] =   N_[iPu242_][i];
                N_Am241_[cellI] =   N_[iAm241_][i];
                N_Am243_[cellI] =   N_[iAm243_][i];
                N_Cm242_[cellI] =   N_[iCm242_][i];
                N_Cm243_[cellI] =   N_[iCm243_][i];
                N_Cm244_[cellI] =   N_[iCm244_][i]; 

                formFactor_[cellI] = Qprofile[i];

                // Calculate Bu as MWd/MTUO2 and MWd/MTU
                // It might duplicate the field Bu, but it can be useful for example 
                // when using "fromList" as "burnupType". In that case the field Bu 
                // represents only the average burnup
                const scalar Buold_ = Bu_av_old[sliceID]*0.8815;
                const scalar Bunew_ = Bu_av[sliceID]*0.8815;
                const scalar BukgUO2_old = Bu_.oldTime().internalField()[cellI];

                scalar newBu = BukgUO2_old + (Bunew_ - Buold_)*formFactor_[cellI];
                Bui[cellI] = newBu;
                BukgU_[cellI] = newBu/0.8815;
            }
        }
    }

    // Bu_.correctBoundaryConditions();
    BukgU_ = Bu_/0.8815;
}    

const Foam::scalarField Foam::burnupLassmann::SigmaA() const
{

    // Compute and return SigmaA_ in (1/m)
    const scalarField SigmaA = 1e2*(
                    sa_th_[iU235_] * N_U235_.internalField()  +
                    sa_th_[iU238_] * N_U238_.internalField()  +
                    sa_th_[iPu239_] * N_Pu239_.internalField() +
                    sa_th_[iPu240_] * N_Pu240_.internalField() +
                    sa_th_[iPu241_] * N_Pu241_.internalField() +
                    sa_th_[iPu242_] * N_Pu242_.internalField() 
                );

    return SigmaA;
}

const Foam::scalarField Foam::burnupLassmann::neutronD() const
{

    // Compute total N in cellI (1/cm3)
    const scalarField Ntot = N_U235_.internalField() +
                             N_U238_.internalField()  +
                             N_Pu239_.internalField() +
                             N_Pu240_.internalField() +
                             N_Pu241_.internalField() +
                             N_Pu242_.internalField() ;

    // Create table of internal radius VS p2 
    
    // scalarField radii("4(0.0e-3, 0.5e-3, 1.0e-3, 1.5e-3)");
    scalarField radiiList({0.0e-3, 0.5e-3, 1.0e-3, 1.5e-3});
    scalarField p2List({0.026759, 0.036073, 0.036558, 0.036084});
    scalarField p3List({-0.497919, -0.433075, -0.423868, -0.421376});

    scalarInterpolateTable p2Table(radiiList, p2List, interpolateTableBase::interpolationMethod::LINEAR);
    scalarInterpolateTable p3Table(radiiList, p3List, interpolateTableBase::interpolationMethod::LINEAR);

    // Build the internal field of D to be returned by this function. D is built
    // material by material
    scalarField p1(mesh_.nCells(), 300e-24);
    scalarField p2(mesh_.nCells(), VSMALL);
    scalarField p3(mesh_.nCells(), VSMALL);
    
    const labelListList& matAddrList(mat_.matAddrList());

    scalarField D(mesh_.nCells(), GREAT);

    forAll(mat_.materialsList(), i)
    {
        const labelList& addr(matAddrList[i]);

        if( isA<fuelMaterial>(mat_.materialsList()[i]) )
        {
            // Initialize local concentration
            forAll(addr, i)
            {
                const label cellI = addr[i];
                // Coefficients for diffusion coefficient
                p2[cellI] =  p2Table(rIn_[cellI]);
                p3[cellI] =  p3Table(rIn_[cellI]);

                // Diffusion coeff. in (m)
                D[cellI] = 1e-2/(3.0*p1[cellI]*Ntot[cellI]
                    *p2[cellI]*pow(rOut_[cellI]-rIn_[cellI], p3[cellI]));
            }
        }  
    }            

    return D;
}

void Foam::burnupLassmann::calcRinRout(const labelList& addr)
{
    // Initialize radii
    scalar ri(GREAT);
    scalar ro(-GREAT);

    // Mesh points
    const pointField& points = mesh_.points();

    forAll(addr, i)
    {
        const label cellI = addr[i];

        forAll( mesh_.cellPoints()[cellI], j )
        {
            const vector pi = points[mesh_.cellPoints()[cellI][j]];
            scalar r = sqrt(pow(pi.x(),2) + pow(pi.y(),2));
            if ( ri > r )
            {
                ri = r;
            }
            if ( ro < r )
            {
                ro = r;
            }  
        }
    }

    // Assign inner and outer radius to the cell
    forAll(addr, i)
    {
        const label cellI = addr[i];
        rIn_[cellI] = ri;
        rOut_[cellI] = ro;

    }
    //Info << "\tInner radius : " << ri  <<  " m" << endl;
    //Info << "\tOuter radius : " << ro <<  " m" << endl;
}

void Foam::burnupLassmann::initNuclearData(const labelList& addr)
{
    // Material enrichment as U235/U (fraction)
    const scalar e = enrichment_[addr[0]];
    scalar e_wt = e*100.0;

    // Clamp to [2,5] so we don't go outside table
    e_wt = max(2.0, min(5.0, e_wt));

    // Lower/upper integer grid points (2,3,4,5)
    const int eLow  = static_cast<int>(floor(e_wt + SMALL));
    const int eHigh = static_cast<int>(ceil(e_wt - SMALL));

    // If exactly at integer, just use that set
    if (eLow == eHigh)
    {
        setXSForIntEnrichment(eLow, sf_, sc_);
    }
    else
    {
        // Linear interpolation weight in [0,1]
        const scalar w = (e_wt - scalar(eLow)) / scalar(eHigh - eLow);
    
        // Temporary storage for low/high XS
        scalarField sfLow(Nelements_, 0.0);
        scalarField scLow(Nelements_, 0.0);
        scalarField sfHigh(Nelements_, 0.0);
        scalarField scHigh(Nelements_, 0.0);
    
        setXSForIntEnrichment(eLow,  sfLow,  scLow);
        setXSForIntEnrichment(eHigh, sfHigh, scHigh);
        
        // Interpolate for all isotopes
        for (int iso = 0; iso < Nelements_; ++iso)
        {
            sf_[iso] = (1.0 - w)*sfLow[iso] + w*sfHigh[iso];
            sc_[iso] = (1.0 - w)*scLow[iso] + w*scHigh[iso];
        }
    }

    // Thermal absorption (used for flux calculation)
    // XS computed at 50 MWd/kgHM, e = 2 wt-%,
    // 2 E-groups w/ default 2 groups division serpent
    sa_th_[iU235_] = 338.43e-24;
    sa_th_[iU238_] = 1.47e-24;
    sa_th_[iPu239_] = 1151.35e-24;
    sa_th_[iPu240_] = 193.85e-24;
    sa_th_[iPu241_] = 1035.21e-24;
    sa_th_[iPu242_] = 11.07e-24;

    // U238(n,2n)Np237 reaction
    // XS computed at 50 MWd/kgHM, e = 2 wt-%,
    // 2 E-groups w/ threshold at 0.35 eV
    sn2n_ = 0.00550912e-24;

}