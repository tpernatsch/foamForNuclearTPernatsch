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

#include "burnupLassmannFBR.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "Tuple2.H"
#include "globalFieldLists.H"
#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(burnupLassmannFBR, 0);
    addToRunTimeSelectionTable
    (
        burnup, 
        burnupLassmannFBR, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::HashTable<Foam::label> Foam::burnupLassmannFBR::setNuclidesIDs()
{
    // Create table of all the isotopes with their index
    // It associates a unique index for each isotope
    Foam::HashTable<Foam::label> table
    {
        {"He4", 0},
        // Uranium isotopes
        {"U234" , 1},
        {"U235" , 2},
        {"U236" , 3},
        {"U237" , 4},
        {"U238" , 5},
        // Neptunium isotopes
        {"Np237", 6},
        {"Np238", 7},
        {"Np239", 8},
        // Plutonium isotopes
        {"Pu238", 9},
        {"Pu239", 10},
        {"Pu240", 11},
        {"Pu241", 12},
        {"Pu242", 13},
        {"Pu243", 14},
        // Americium isotopes
        {"Am241", 15},
        {"Am242", 16},
        {"Am242m", 17},
        {"Am243", 18},
        {"Am244", 19},
        // Curium isotopes
        {"Cm242", 20},
        {"Cm243", 21},
        {"Cm244", 22},
        {"Cm245", 23}
    };

    return table;
}

Foam::HashTable<Foam::scalar> Foam::burnupLassmannFBR::setNuclidesLambdas()
{
    // List of half lives in years
    Foam::HashTable<Foam::scalar> table
    {
        // He4
        {"He4",    0},
        // Uranium isotopes
        {"U234" ,  0},
        {"U235" ,  0},
        {"U236" ,  0},
        {"U237" ,  log(2.) / (6.750)},
        {"U238" ,  0},
        // Neptunium isotopes
        {"Np237",  0},
        {"Np238",  log(2.) / (2.117)},
        {"Np239",  log(2.) / (2.36)},
        // Plutonium isotopes
        {"Pu238",  log(2.) / (8.77e+01  * 365.)},
        {"Pu239",  0},
        {"Pu240",  0},
        {"Pu241",  log(2.) / (1.43e+01  * 365.)},
        {"Pu242",  log(2.) / (3.73e+05  * 365.)},
        {"Pu243",  log(2.) / (4.956     / 24.)}, // ! Very short lived
        // Americium isotopes
        {"Am241",  log(2.) / (432.2   * 365.)}, 
        {"Am242",  log(2.) / (16.02   / 24.)}, // ! Very short lived
        {"Am242m", log(2.) / (141.0   * 365.)},
        {"Am243",  log(2.) / (7350.0  * 365.)},
        {"Am244",  log(2.) / (10.0    / 24.)}, // ! Very short lived
        // Curium isotopes
        {"Cm242",  log(2.) / (162.8  * 1.00)},
        {"Cm243",  log(2.) / ( 29.1  * 365.)},
        {"Cm244",  log(2.) / ( 18.1  * 365.)},
        {"Cm245",  log(2.) / (8.5e3  * 365.)}
    };

    return table;
}

void Foam::burnupLassmannFBR::createNuclidesFields()
{
    // Initialize hashtable elements names (sum of all isotopes)
    HashTable<label> elemNames; 
    

    forAllIter(HashTable<label>, ID_, iter)
    {
        // Get name of this isotope 
        word name(iter.key());

        // Assign isotope number density field
        nuclPtrList_.set
        (
            iter(),
            &createOrLookup<scalar>
            (
                mesh_,
                Foam::word("N_")+name,
                dimless/dimVolume,
                0.0,
                "zeroGradient"
            )
        );

        // Get element name of this isotope removing num characters, e.g. Pu239 -> Pu
        std::string elemName = name;
        // Erase everything from first digit onwards
        elemName.erase
        (
            std::find_if(elemName.begin(), elemName.end(), ::isdigit), 
            elemName.end()
        );

        // If element not found in elemNames, append it
        if (not elemNames.found("N_"+elemName))
        {
            elemNames.insert("N_"+elemName, 1); // 1 is a dummy value
        }
    }

    // Create all elements field in elemNames hasTable
    forAllIter(HashTable<label>, elemNames, iter)
    {
        elemPtrList_.append
        (
            &createOrLookup<scalar>
            (
                mesh_,
                iter.key(),
                dimless/dimVolume,
                0.0,
                "zeroGradient"
            )
        );
    }   
}

void Foam::burnupLassmannFBR::setXS()
{
    // Get type of reactor
    const word type = globOpt_.reactorType();

    // Initialize quantity with which XS are parametrized, it can be U-235 
    // enrichment or Pu content, depending on reactorType
    scalarField parameter;

    // Initialize quantity with which the XS tables are later interrogated, it 
    // can be U-235 enrichment or Pu content, depending on reactorType
    scalar x(-1);

    if ( type == "LWR" )
    {
        // For LWRs, XS are expected to be provided in constant/XSdict, 
        // parametrized as a function of the U-235 enrichment
        parameter = scalarField(XSdict_.lookup("enrichments"));

        // 
        // NOTICE, this burnup class burns only one material, 'cause only one
        // set of XS is now supported. To be changed 
        // 
        forAll(mat_.materialsList(), i)
        {
            bool isFuel(isA<fuelMaterial>(mat_.materialsList()[i]));

            if(isFuel)
            {
                const UO2& UO2Mat = 
                    refCast<const UO2>(mat_.materialsList()[i]);
                
                x = UO2Mat.enrichment();
            }
            break;
        }

        // Set mean energy per fission event (in MWd)
        Efiss_ = 3.61e-22;

        // Set Branching ratios
        BR_ng_Am242_ = 0.85;
        BR_Am241_ = 0.919;   // src: Serpent Wiki page
        BR_Am243_ = 0.06260; // src: Serpent Wiki page
    }
    else if ( type == "FBR" )
    {
        // For FBRs, XS are expected to be provided in constant/XSdict, 
        // parametrized as a function of the Pu wt.% content
        parameter = scalarField(XSdict_.lookup("PuValues"));

        // 
        // NOTICE, this burnup class burns only one material, 'cause only one
        // set of XS is now supported. To be changed 
        // 
        forAll(mat_.materialsList(), i)
        {
            bool isFuel(isA<fuelMaterial>(mat_.materialsList()[i]));

            if(isFuel)
            {
                const UPuO2& UPuO2Mat = 
                    refCast<const UPuO2>(mat_.materialsList()[i]);
                
                x = UPuO2Mat.PuMetalRatio();
            }
            break;
        }

        // Set mean energy per fission event (in MWd)
        Efiss_ = 3.80e-22;

        // Set Branching ratios
        BR_ng_Am242_ = 0.88;
        // BR_Am241_ = 0.8419; // src: "Hiruta et Al., Impact of Americium-241..."
        BR_Am241_ = 0.919; // src: Serpent Wiki page
        BR_Am243_ = 0.06260; // src: Serpent Wiki page
    }
    else
    {
        x = -1;

        FatalErrorInFunction()
        << "Unknown reactor type " << nl
        << "Valid types are:" << endl
        << "  - LWR " << nl 
        << "  - FBR " << nl 
        << exit(FatalError);
    }

    // Read fission XS data as field of fields
    FieldField<Field, scalar> fissionXSdata = 
        PtrList<scalarField>
        (
            XSdict_.lookup("fissionXS"), scalarFieldFieldINew()
        );

    // Read capture XS data as field of fields
    FieldField<Field, scalar> captureXSdata = 
        PtrList<scalarField>
        (
            XSdict_.lookup("captureXS"), scalarFieldFieldINew()
        );

    // Read (n,2n) U238 XS data as scalarField
    scalarField n2nXSdata = XSdict_.lookup("n2nU238");

    // Read (n,3n) U238 XS data as scalarField
    scalarField n3nXSdata = XSdict_.lookup("n3nU238");

    // Read (n,alpha) O16 XS data as scalarField
    scalarField nAlphaXSdata = XSdict_.lookup("nAlphaO16");

    // Check that data are consistent
    forAll(fissionXSdata, i)
    {
        if
        (
            fissionXSdata[i].size() != ID_.size()
            or
            captureXSdata[i].size() != ID_.size()
        )
        {
            FatalErrorInFunction()
            << "XS table dimension inconsistent with number of nuclides. " << nl
            << "There should be " << ID_.size() << " entries per XS list."
            << "XS are needed for the following isotopes : " << nl
            << ID_.sortedToc() << nl 
            << "NOTE : XS must be provided in the table in order of ascending " 
            << "atomic number Z ! I.e. He, U, Np, Pu, Am etc.." << nl
            << exit(FatalError);
        }
    }
    
    // Create fission XS interpolation table
    scalarFieldInterpolateTable fissionXStable
        (
            parameter,
            fissionXSdata,
            interpolateTableBase::interpolationMethodNames_["linear"],
            interpolateTableBase::outOfBoundsMethodNames_["fixed"]
        );
    
    // Create capture XS interpolation table
    scalarFieldInterpolateTable captureXStable
        (
            parameter,
            captureXSdata,
            interpolateTableBase::interpolationMethodNames_["linear"],
            interpolateTableBase::outOfBoundsMethodNames_["fixed"]
        );

    // Create (n,2n) U238 XS interpolation table
    scalarInterpolateTable n2nXStable
        (
            parameter,
            n2nXSdata,
            interpolateTableBase::interpolationMethodNames_["linear"],
            interpolateTableBase::outOfBoundsMethodNames_["fixed"]
        );


    // Create (n,3n) U238 XS interpolation table
    scalarInterpolateTable n3nXStable
        (
            parameter,
            n3nXSdata,
            interpolateTableBase::interpolationMethodNames_["linear"],
            interpolateTableBase::outOfBoundsMethodNames_["fixed"]
        );

    // Create (n,alpha) O16 XS interpolation table
    scalarInterpolateTable nAlphaXStable
        (
            parameter,
            nAlphaXSdata,
            interpolateTableBase::interpolationMethodNames_["linear"],
            interpolateTableBase::outOfBoundsMethodNames_["fixed"]
        );

    // Assign fission and capture XS (convert to barns)
    sf_ = fissionXStable(x) * 1e-24;
    sc_ = captureXStable(x) * 1e-24;
    sa_ = sf_ + sc_;
    sn2nU238_ = n2nXStable(x) * 1e-24;
    sn3nU238_ = n3nXStable(x) * 1e-24;
    snAlphaO16_ = nAlphaXStable(x) * 1e-24;

    // Assign thermal absorption XS only in LWR case, otherwise leave to zero
    if ( type == "LWR" )
    {
        // Thermal absorption (used for flux calculation)
        // XS computed at 50 MWd/kgHM, e = 2 wt-%,
        // 2 E-groups w/ default 2 groups division serpent
        sa_th_[ID_["U235"]] = 338.43e-24;
        sa_th_[ID_["U238"]] = 1.47e-24;
        sa_th_[ID_["Pu239"]] = 1151.35e-24;
        sa_th_[ID_["Pu240"]] = 193.85e-24;
        sa_th_[ID_["Pu241"]] = 1035.21e-24;
        sa_th_[ID_["Pu242"]] = 11.07e-24;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burnupLassmannFBR::burnupLassmannFBR
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& burnupDict
)
:
    burnupFromPower(mesh, mat, burnupDict),
    globOpt_(mesh_.lookupObject<globalOptions>("globalOptions")),
    XSdict_
    (
        IOobject
        (
            "XSdata",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
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
        dimensionedScalar("formFactor", dimless, 1),    
        "zeroGradient"
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
        dimensionedScalar("Bu_kgU", dimless, 0),
        "zeroGradient"
    ),
    ID_(setNuclidesIDs()),
    nuclLambdas_(setNuclidesLambdas()),
    nuclPtrList_(ID_.size()),
    elemPtrList_(),
    sf_(ID_.size(),0.0),
    sc_(ID_.size(),0.0),
    sa_(ID_.size(),0.0),
    sa_th_(ID_.size(),0.0),
    sn2nU238_(0.0),
    sn3nU238_(0.0),
    BR_beta_Am242_(0.827),
    BR_ng_Am242_(-GREAT),
    BR_Am241_(-GREAT),
    BR_Am243_(-GREAT),
    snAlphaO16_(0.0),
    HMoverFuelRatio_(0.8815),
    Efiss_(-GREAT),
    rIn_(mesh_.nCells(),GREAT),
    rOut_(mesh_.nCells(),-GREAT),
    mapper_(nullptr)
{   
    // Fill pointer list with nuclides volScalarFields
    createNuclidesFields();

    // Set cross section fields (read from dictionary + interpolate)
    setXS();

    // Read precision of convergence loop for Bateman's equations
    precision_ = burnupDict_.lookupOrDefault("convergencePrecision", 1e-2);

    // Compute inner and outer radii (for radial profile resonance capture)
    const labelListList& matAddrList(mat_.matAddrList());
    forAll(mat_.materialsList(), i)
    {   
        const labelList& addr(matAddrList[i]);
        calcRinRout(addr);
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::burnupLassmannFBR::~burnupLassmannFBR()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::burnupLassmannFBR::correct()
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
    
    // Get slice addr list
    const PtrList<labelList>& sliceAddr(mapper_->sliceAddrList());

    // Ref to power density field
    const volScalarField& Q = mesh_.lookupObject<volScalarField>("Q");

    // Compute Q volume-average on each slice
    const scalarField& QAvg(mapper_->sliceAverage(Q));

    // Ref to density field
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    // Compute rho volume-average on each slice
    const scalarField& rhoAvg(mapper_->sliceAverage(rho));

    // Ref to O16 concentration (created by fuelMaterial)
    const volScalarField& N_O16 = mesh_.lookupObject<volScalarField>("N_O16");

    // Compute Bu volume-average on each slice and convert
    // from MWd/ton -> MWd/tonHM
    scalarField BuAvg(mapper_->sliceAverage(Bu_));
    scalarField BuAvgOld(mapper_->sliceAverage(Bu_.oldTime()));
    BuAvg /= HMoverFuelRatio_;
    BuAvgOld /= HMoverFuelRatio_;

    // Reference to volumes
    const scalarField& V(mesh_.V());

    // Compute vol per each slice
    scalarField volPerSlice(nSlices, 0.0);
    forAll(sliceAddr, sliceID)
    {  
        const labelList& addr(sliceAddr[sliceID]);
        scalarField Vslice(V, addr);
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
           
            dBu[cellI] = (BuAvg[sliceID] - BuAvgOld[sliceID]) / 1000;
        }
    }

// READ INPUT AND OLD TIME CONCENTRATIONS 

    const List<bool>& isFuelList(mapper_->isFuel());
    forAll(sliceAddr, sliceID)
    {
        if(isFuelList[sliceID])
        {
            const List<label>& addrSlice_(sliceAddr[sliceID]);
            scalar localSize(addrSlice_.size());
            scalarField Vlocal(V, addrSlice_);

            // Get average power density for this slice in MW/m3
            const scalar QAvgI = QAvg[sliceID] * 1e-6;

            // Get average density for this slice
            const scalar rhoAvgI = rhoAvg[sliceID];

            // Compute average HM density
            const scalar rhoHM = rhoAvgI * HMoverFuelRatio_;

            // Init old time concentrations [1/cm3]
            scalarRectangularMatrix Nold(ID_.size(),localSize);
            scalarField NavgOld(ID_.size(),0.0);
            scalarField Navg(ID_.size(),0.0);
            forAll(addrSlice_, i)
            {
                const label cellI = addrSlice_[i];

                forAllIter(HashTable<label>, ID_, iter)
                {
                    Nold[iter()][i]  = nuclPtrList_[iter()].oldTime()()[cellI];
                    NavgOld[iter()] += Nold[iter()][i]*Vlocal[i];
                }
            }
            reduce(NavgOld, sumOp<scalarField>());
            NavgOld /= volPerSlice[sliceID];
            Navg = NavgOld;

    // NUCLEAR DATA

            // Branching ratios
            const scalar BRPu2 = (1-BR_beta_Am242_) * BR_ng_Am242_; 
            const scalar BRCm2 = BR_beta_Am242_ * BR_ng_Am242_;
            
            // Ternary fission yield (from paper Cechet et al.)
            const scalar ternFissY = 0.22e-2;

            // Fission rate density (fiss./cm3/day)
            const scalar Fdot = (QAvgI*1e-6)/Efiss_;

    // RADIAL SHAPE FUNCTIONS FOR RESONANCE CAPTURE

            scalarField frU(localSize, 0.0);
            scalarField frPu(localSize, 0.0);

            if ( globOpt_.reactorType() == "LWR" )
            {
                scalar normU = 0;
                scalar normPu = 0;

                const volVectorField& C(mesh_.C());

                // Calculate volume-averaged radial shape functions
                forAll(addrSlice_, i)
                {
                    const label cellI = addrSlice_[i];

                    scalar r = pow(pow(C[cellI][0],2)+pow(C[cellI][1],2),0.5);

                    frU[i]  = (1 + 3.45*exp(-3.0*pow((rOut_[cellI] - r)*1e3, 0.45 )));
                    frPu[i] = (1 + 0.94*exp(-3.0*pow((rOut_[cellI] - r)*1e3, 0.45 )));

                    // Calculate average frU and frPu
                    normU += frU[i]*Vlocal[i];
                    normPu += frPu[i]*Vlocal[i];
                }

                reduce(normU, sumOp<scalar>());
                reduce(normPu, sumOp<scalar>());

                normU /= volPerSlice[sliceID];
                normPu /= volPerSlice[sliceID];

                frU = frU/normU;
                frPu = frPu/normPu;

            }
            else if ( globOpt_.reactorType() == "FBR" )
            {
                // Neglect resonance capture
                frU = scalarField(localSize, 1.0);
                frPu = scalarField(localSize, 1.0);
            }
            else
            {
                FatalErrorInFunction()
                << "Unknown reactor type " << nl
                << "Valid types are:" << endl
                << "  - LWR " << nl 
                << "  - FBR " << nl 
                << exit(FatalError);
            }

    // CONVERGENCE LOOP      

            scalarRectangularMatrix NoldIter(ID_.size(), localSize);
            scalarRectangularMatrix N(ID_.size(), localSize);
            scalarField Qprofile(localSize, 0.0);  
            N = Nold;

            scalarField dBuLocal(dBu, addrSlice_);
            scalarField dBuLocalNoProf(dBu, addrSlice_);

            // Compute macro fission XS
            scalar SFtot = 0;
            scalar SFtotOld = 0;
            scalarField SF(localSize,0.0);
            forAllIter(HashTable<label>, ID_, iter)
            {
                SFtot += sf_[iter()] * Navg[iter()];
                SFtotOld += sf_[iter()] * NavgOld[iter()];

                forAll(addrSlice_, i)
                {
                    SF[i] += sf_[iter()]*N[iter()][i];
                }
            }

            SFtot = SFtotOld;

            // Get delta T for radioactive decay in days
            const scalar& deltaT = mesh_.time().deltaTValue() / 3600 / 24;
            
            scalar difference = GREAT;
            while (difference > precision_)
            {
                NoldIter = N;
                
                scalarField S(ID_.size(),0.0);

                // Matrix creation and new compositions calculation
                forAll(addrSlice_, i)
                {
                    // Compute A from Eq.3 of Lassmann paper
                    scalar A = 1e-6 * rhoHM / ( Efiss_*SF[i] );

                    // Init matrix of coefficients for Bateman's eq.
                    scalarSquareMatrix M(ID_.size(),0.0);
                    // Init array of unknowns for Bateman's eq.
                    scalarList X(ID_.size(),0.0);
                    // Init array of old-time Ns for Bateman's eq.
                    scalarList B(ID_.size(),0.0);
                    // Init identity matrix
                    scalarSquareMatrix I(ID_.size(),0.0);

                    // Init matrix for absorption contribution
                    scalarSquareMatrix ABS(ID_.size(),0.0);
                    // Init matrix for capture contribution
                    scalarSquareMatrix CAP(ID_.size(),0.0);
                    // Init matrices for decay source and sink
                    scalarSquareMatrix DEC_src(ID_.size(),0.0);
                    scalarSquareMatrix DEC_sink(ID_.size(),0.0);

                    forAllIter(HashTable<label>, ID_, iter)
                    {
                        const word name = iter.key();

                        // Set identity matrix diagonal
                        I[iter()][iter()] = 1.0;

                        // Set decay sink (this is a diagonal term)
                        DEC_sink[iter()][iter()] = nuclLambdas_[name] * deltaT;

                        // Set old time concentration
                        B[iter()] = Nold[iter()][i];

                        // Set absorption contribution (this is a diagonal term)
                        if (name == "U238")
                        {
                            ABS[iter()][iter()]=
                            ( frU[i]*sa_[iter()] + sn2nU238_ ) * A * dBuLocal[i];
                        }
                        else if (name == "Pu240")
                        {
                            ABS[iter()][iter()]=
                                frPu[i] * sa_[iter()] * A * dBuLocal[i];
                        }
                        else
                        {
                            ABS[iter()][iter()] = sa_[iter()]*A*dBuLocal[i];
                        }
                    }

                    // Set capture contribution (non-diagonal)
                    CAP[ ID_["U235"]  ][ ID_["U234"]  ] = -1 * sc_[ID_["U234"]]  * A * dBuLocal[i];
                    CAP[ ID_["U236"]  ][ ID_["U235"]  ] = -1 * sc_[ID_["U235"]]  * A * dBuLocal[i];
                    CAP[ ID_["U237"]  ][ ID_["U236"]  ] = -1 * sc_[ID_["U236"]]  * A * dBuLocal[i];
                    CAP[ ID_["U238"]  ][ ID_["U237"]  ] = -1 * sc_[ID_["U237"]]  * A * dBuLocal[i];
                   
                    CAP[ ID_["Np237"] ][ ID_["U236"]  ] = -1  * sc_[ID_["U236"]]  * A * dBuLocal[i];
                    CAP[ ID_["Np238"] ][ ID_["Np237"] ] = -1 * sc_[ID_["Np237"]]  * A * dBuLocal[i];
                    CAP[ ID_["Np239"] ][ ID_["U238"]  ] = -1  * sc_[ID_["U238"]]  * A * dBuLocal[i];
                    CAP[ ID_["Np239"] ][ ID_["Np238"] ] = -1 * sc_[ID_["Np238"]]  * A * dBuLocal[i];
                   
                    CAP[ ID_["Pu239"] ][ ID_["Pu238"] ] = -1       * sc_[ID_["Pu238"]] * A * dBuLocal[i];
                    CAP[ ID_["Pu239"] ][ ID_["U238"]  ] = -frU[i]  * sc_[ID_["U238"]]  * A * dBuLocal[i];
                    CAP[ ID_["Pu240"] ][ ID_["Pu239"] ] = -1       * sc_[ID_["Pu239"]] * A * dBuLocal[i];
                    CAP[ ID_["Pu241"] ][ ID_["Pu240"] ] = -frPu[i] * sc_[ID_["Pu240"]] * A * dBuLocal[i];
                    CAP[ ID_["Pu242"] ][ ID_["Pu241"] ] = -1       * sc_[ID_["Pu241"]] * A * dBuLocal[i];
                    CAP[ ID_["Pu243"] ][ ID_["Pu242"] ] = -1       * sc_[ID_["Pu242"]] * A * dBuLocal[i];
                    
                    CAP[ ID_["Am242"] ][ ID_["Am241"] ] = -BR_Am241_ * sc_[ID_["Am241"]] * A * dBuLocal[i];
                    CAP[ ID_["Am242m"] ][ ID_["Am241"] ] = -(1-BR_Am241_) * sc_[ID_["Am241"]] * A * dBuLocal[i];
                    CAP[ ID_["Am243"] ][ ID_["Am242"] ] = - 1 * sc_[ID_["Am242"]] * A * dBuLocal[i];
                    CAP[ ID_["Am243"] ][ ID_["Am242m"] ] = - 1 * sc_[ID_["Am242m"]] * A * dBuLocal[i];
                    CAP[ ID_["Am244"] ][ ID_["Am243"] ] = -BR_Am243_ * sc_[ID_["Am243"]] * A * dBuLocal[i];
                    
                    CAP[ ID_["Cm243"] ][ ID_["Cm242"] ] = -1 * sc_[ID_["Cm242"]] * A * dBuLocal[i];
                    CAP[ ID_["Cm244"] ][ ID_["Cm243"] ] = -1 * sc_[ID_["Cm243"]] * A * dBuLocal[i];
                    CAP[ ID_["Cm245"] ][ ID_["Am244"] ] = -1 * sc_[ID_["Am244"]] * A * dBuLocal[i];
                    CAP[ ID_["Cm245"] ][ ID_["Cm244"] ] = -1 * sc_[ID_["Cm244"]] * A * dBuLocal[i];

                    // Not a capture: U238 + n -> U237 + 2n
                    CAP[ ID_["U237"] ][ ID_["U238"]  ] = -1 * sn2nU238_ * A * dBuLocal[i];

                    // Not a capture: U238 + n -> U236 + 3n
                    CAP[ ID_["U236"] ][ ID_["U238"]  ] = -1 * sn3nU238_ * A * dBuLocal[i];

                    // Set decay source contribution (non-diagonal)
                    DEC_src[ID_["U234"]][ID_["Pu238"]] = -nuclLambdas_["Pu238"] * deltaT;
                    DEC_src[ID_["U236"]][ID_["Pu240"]] = -nuclLambdas_["Pu240"] * deltaT;
                    DEC_src[ID_["Np237"]][ID_["Am241"]] = -nuclLambdas_["Am241"] * deltaT;
                    DEC_src[ID_["Np237"]][ID_["U237"]] = -nuclLambdas_["U237"] * deltaT;
                    DEC_src[ID_["Np239"]][ID_["Am243"]] = -nuclLambdas_["Am243"] * deltaT;
                    DEC_src[ID_["Pu238"]][ID_["Np238"]] = -nuclLambdas_["Np238"] * deltaT;
                    DEC_src[ID_["Pu238"]][ID_["Cm242"]] = -nuclLambdas_["Cm242"] * deltaT;
                    DEC_src[ID_["Pu240"]][ID_["Cm244"]] = -nuclLambdas_["Cm244"] * deltaT;
                    DEC_src[ID_["Pu241"]][ID_["Cm245"]] = -nuclLambdas_["Cm245"] * deltaT;
                    DEC_src[ID_["Pu242"]][ID_["Am242"]] = -(1-BR_beta_Am242_) * nuclLambdas_["Am242"] * deltaT;
                    DEC_src[ID_["Am241"]][ID_["Pu241"]] = -nuclLambdas_["Pu241"] * deltaT;
                    DEC_src[ID_["Am243"]][ID_["Pu243"]] = -nuclLambdas_["Pu243"] * deltaT;
                    DEC_src[ID_["Am242"]][ID_["Am242m"]] = -nuclLambdas_["Am242m"] * deltaT;
                    DEC_src[ID_["Cm242"]][ID_["Am242"]] = -BR_beta_Am242_ * nuclLambdas_["Am242"] * deltaT;
                    
                    // Important: More than 90% of the captures of Am243 generate
                    // Am244m, which has a short half life (~26 min)
                    // Since I don't care abouot tracking the nuclides of Am244m
                    // I just weight the half lives with branching ratios
                    const scalar lambdaAm244metastable = log(2.0)/(26./60./24.);
                    DEC_src[ID_["Cm244"]][ID_["Am244"]] = - 1.0 * 
                    (
                        BR_Am243_ * nuclLambdas_["Am244"]
                        + 
                        (1-BR_Am243_) * lambdaAm244metastable
                    )* deltaT;

                    // He4 production by alpha decay
                    DEC_src[ID_["He4"]][ID_["U234"]] = -nuclLambdas_["U234"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["U235"]] = -nuclLambdas_["U235"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["U236"]] = -nuclLambdas_["U236"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["U238"]] = -nuclLambdas_["U238"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Np237"]] = -nuclLambdas_["Np237"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Pu238"]] = -nuclLambdas_["Pu238"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Pu239"]] = -nuclLambdas_["Pu239"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Pu240"]] = -nuclLambdas_["Pu240"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Pu242"]] = -nuclLambdas_["Pu242"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Am241"]] = -nuclLambdas_["Am241"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Am243"]] = -nuclLambdas_["Am243"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Cm242"]] = -nuclLambdas_["Cm242"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Cm243"]] = -nuclLambdas_["Cm243"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Cm244"]] = -nuclLambdas_["Cm244"] * deltaT;
                    DEC_src[ID_["He4"]][ID_["Cm245"]] = -nuclLambdas_["Cm245"] * deltaT;

                    // He4 production by ternary fission (in source vector)
                    B[ID_["He4"]] += ternFissY * Fdot * deltaT;

                    // He4 production by O16 + n -> C13 + He4 (in source vector)
                    B[ID_["He4"]] += snAlphaO16_ * N_O16[i] * A * dBuLocal[i];

                    // Build matrix of coefficients for Bateman's equations 
                    M = I + ABS + CAP + DEC_sink + DEC_src;
                    
                    // Solve system of equations
                    solve(X, M, B);

                    SF[i] = 0;
                    forAllIter(HashTable<label>, ID_, iter)
                    {
                        N[iter()][i] = X[iter()];
                        S[iter()] += N[iter()][i]*Vlocal[i];
                        SF[i] += sf_[iter()]*N[iter()][i];
                    }
                }

                reduce(S, sumOp<scalarField>());
                S /= volPerSlice[sliceID];

                Navg = max(S, SMALL);
                
                // Flux, burnup and power profiles calculation
                SFtot = 0;
                forAllIter(HashTable<label>, ID_, iter)
                {
                    SFtot += sf_[iter()] * Navg[iter()];
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
                    FatalErrorIn("Foam::burnupLassmannFBR::correct()")
                    << "This burnup model needs a neutronics subsolver to "
                    << "be set. " << exit(FatalError);
                }

                // Compute power profile and normalize it by its average
                Qprofile = fluxProfile * SF;
                Qprofile /= gSum(Qprofile * Vlocal)/volPerSlice[sliceID];

                // Apply Q profile to local bu
                dBuLocal = Qprofile*dBuLocalNoProf;

                // Difference calculation for convergence  
                scalar dif(0);
                scalar maxdif(0);
                forAllIter(HashTable<label>, ID_, iter)
                {
                    forAll(addrSlice_, i)
                    {
                        // const label cellI = addrSlice_[i];
                        
                        if (NoldIter[iter()][i]>SMALL)
                        {
                            dif=mag((N[iter()][i]-NoldIter[iter()][i])/NoldIter[iter()][i]);

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

                // Write nuclides density fields
                forAllIter(HashTable<label>, ID_, iter)
                {
                    nuclPtrList_[iter()][cellI] = N[iter()][i];
                }
                
                formFactor_[cellI] = Qprofile[i];

                // Calculate Bu as MWd/MTUO2 and MWd/MTU
                // It might duplicate the field Bu, but it can be useful for example 
                // when using "fromList" as "burnupType". In that case the field Bu 
                // represents only the average burnup
                const scalar Buold = BuAvgOld[sliceID]*HMoverFuelRatio_;
                const scalar Bunew = BuAvg[sliceID]*HMoverFuelRatio_;
                const scalar BukgUO2_old = Bu_.oldTime().internalField()[cellI];

                scalar newBu = BukgUO2_old + (Bunew - Buold)*formFactor_[cellI];
                Bu_.ref()[cellI] = newBu;
                BukgU_[cellI] = newBu/HMoverFuelRatio_;
            }

        }
    }
    
    // Write elements (sum of all nuclides belonging to that element)
    forAll(elemPtrList_, i)
    {
        const word elemName = elemPtrList_[i].name();
        
        elemPtrList_[i] *= 0.0;

        forAll(nuclPtrList_, j)
        {
            const word nuclName = nuclPtrList_[j].name();

            // Check if elemName is contained in nuclName
            if (nuclName.find(elemName) != std::string::npos) 
            {
                elemPtrList_[i] += nuclPtrList_[j];
            }
        }
    }

    // Bu_.correctBoundaryConditions();
    BukgU_ = Bu_/HMoverFuelRatio_;
}    

const Foam::scalarField Foam::burnupLassmannFBR::SigmaA() const
{
    // Compute and return SigmaA_ in (1/m)
    const scalarField SigmaA = 1e2*(
    sa_th_[ID_["U235"]] * nuclPtrList_[ID_["U235"]]()  +
    sa_th_[ID_["U238"]] * nuclPtrList_[ID_["U238"]]()  +
    sa_th_[ID_["Pu239"]] * nuclPtrList_[ID_["Pu239"]]() +
    sa_th_[ID_["Pu240"]] * nuclPtrList_[ID_["Pu240"]]() +
    sa_th_[ID_["Pu241"]] * nuclPtrList_[ID_["Pu241"]]() +
    sa_th_[ID_["Pu242"]] * nuclPtrList_[ID_["Pu242"]]() 
    );

    return SigmaA;
}

const Foam::scalarField Foam::burnupLassmannFBR::neutronD() const
{

    // Compute total N in cellI (1/cm3)
    const scalarField Ntot = nuclPtrList_[ID_["U235"]]() +
                             nuclPtrList_[ID_["U238"]]()  +
                             nuclPtrList_[ID_["Pu239"]]() +
                             nuclPtrList_[ID_["Pu240"]]() +
                             nuclPtrList_[ID_["Pu241"]]() +
                             nuclPtrList_[ID_["Pu242"]]();

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

void Foam::burnupLassmannFBR::calcRinRout(const labelList& addr)
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
