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

#include "fgrSCIANTIX.H"
#include "addToRunTimeSelectionTable.H"
#include "gapGasModel.H"
#include "zeroGradientFvPatchFields.H"
#include "globalFieldLists.H"
#include <string>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fgrSCIANTIX, 0);
    addToRunTimeSelectionTable(fissionGasRelease, fgrSCIANTIX, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fgrSCIANTIX::readSCIANTIXinput()
{
    // Construct the file path
    Foam::fileName inputCheckFilePath = mesh_.time().path() / "input_check.txt";

    // Create output file for checking input reading
    OFstream inputCheckFile(inputCheckFilePath);

    // Assign SCIANTIX inputs to sciantix options
    Sciantix_options[0] = inputSCIANTIX_.lookupOrDefault("iverification",0);
    Sciantix_options[1] = inputSCIANTIX_.lookupOrDefault("igrain_growth",0);
    Sciantix_options[2] = inputSCIANTIX_.lookupOrDefault("iinert_gas_behavior",0);
    Sciantix_options[3] = inputSCIANTIX_.lookupOrDefault("igas_diffusion_coefficient",0);
    Sciantix_options[4] = inputSCIANTIX_.lookupOrDefault("iintra_bubble_evolution",0);
    Sciantix_options[5] = inputSCIANTIX_.lookupOrDefault("ibubble_radius",0);
    Sciantix_options[6] = inputSCIANTIX_.lookupOrDefault("iresolution_rate",0);
    Sciantix_options[7] = inputSCIANTIX_.lookupOrDefault("itrapping_rate",0);
    Sciantix_options[8] = inputSCIANTIX_.lookupOrDefault("inucleation_rate",0);
    Sciantix_options[9] = inputSCIANTIX_.lookupOrDefault("isolver",0);
    Sciantix_options[10] = inputSCIANTIX_.lookupOrDefault("iformat_output",0);
    Sciantix_options[11] = inputSCIANTIX_.lookupOrDefault("igrain_boundary_vacancy_diffusion_coefficient",0);
    Sciantix_options[12] = inputSCIANTIX_.lookupOrDefault("igrain_boundary_behaviour",0);
    Sciantix_options[13] = inputSCIANTIX_.lookupOrDefault("igrain_boundary_micro_cracking",0);
    Sciantix_options[14] = inputSCIANTIX_.lookupOrDefault("igrain_recrystallization",0);
    Sciantix_options[15] = inputSCIANTIX_.lookupOrDefault("ifuel_reactor_type",0);
    Sciantix_options[16] = inputSCIANTIX_.lookupOrDefault("igas_effective_coefficient",0);
    Sciantix_options[17] = inputSCIANTIX_.lookupOrDefault("igas_sweeping",0);
    Sciantix_options[18] = inputSCIANTIX_.lookupOrDefault("imicro_cracking_span",0);
    // Input scaling factors
    Sciantix_scaling_factors[0] = inputSCIANTIX_.lookupOrDefault("sf_resolution_rate", 1.0);
    Sciantix_scaling_factors[1] = inputSCIANTIX_.lookupOrDefault("sf_trapping_rate", 1.0);
    Sciantix_scaling_factors[2] = inputSCIANTIX_.lookupOrDefault("sf_nucleation_rate", 1.0);
    Sciantix_scaling_factors[3] = inputSCIANTIX_.lookupOrDefault("sf_diffusion_rate", 1.0);

    // Write to check file
    inputCheckFile << "iverification = " << Sciantix_options[0] << endl;
    inputCheckFile << "igrain_growth = " << Sciantix_options[1] << endl;
    inputCheckFile << "iinert_gas_behavior = " << Sciantix_options[2] << endl;
    inputCheckFile << "igas_diffusion_coefficient = " << Sciantix_options[3] << endl;
    inputCheckFile << "iintra_bubble_evolution = " << Sciantix_options[4] << endl;
    inputCheckFile << "ibubble_radius = " << Sciantix_options[5] << endl;
    inputCheckFile << "iresolution_rate = " << Sciantix_options[6] << endl;
    inputCheckFile << "itrapping_rate = " << Sciantix_options[7] << endl;
    inputCheckFile << "inucleation_rate = " << Sciantix_options[8] << endl;
    inputCheckFile << "isolver = " << Sciantix_options[9] << endl;
    inputCheckFile << "iformat_output = " << Sciantix_options[10] << endl;
    inputCheckFile << "igrain_boundary_vacancy_diffusion_coefficient = " << Sciantix_options[11] << endl;
    inputCheckFile << "igrain_boundary_behaviour = " << Sciantix_options[12] << endl;
    inputCheckFile << "igrain_boundary_micro_cracking = " << Sciantix_options[13] << endl;
    inputCheckFile << "igrain_recrystallization = " << Sciantix_options[14] << endl;
    inputCheckFile << "ifuel_reactor_type = " << Sciantix_options[15] << endl;
    inputCheckFile << "igas_effective_coefficient = " << Sciantix_options[16] << endl;
    inputCheckFile << "igas_sweeping = " << Sciantix_options[17] << endl;
    inputCheckFile << "imicro_cracking_span = " << Sciantix_options[18] << endl;

    for(int n = 0; n < 2; n++)
    {
        inputCheckFile << Time_input[n] << "\t";
        inputCheckFile << Temperature_input[n] << "\t";
        inputCheckFile << Fissionrate_input[n] << "\t";
        inputCheckFile << Hydrostaticstress_input[n] << "\t";
        inputCheckFile << endl;
    }

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fgrSCIANTIX::fgrSCIANTIX
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
    refabricatedRod_
    (       
#ifdef OPENFOAMFOUNDATION
        fissionGasReleaseDict_.isNull() ? 
#elif OPENFOAMESI
        fissionGasReleaseDict_.isNullDict() ? 
#endif             
        false :
        fissionGasReleaseDict_.lookupOrDefault("refabricatedRod", false)
    ),
    inputSCIANTIX_(fissionGasReleaseDict_.subOrEmptyDict("SCIANTIX")),
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
        (mesh, "SCIANTIX::Gas_released_t0", sciantixFieldsDict_, 
            (refabricatedRod_ & addToRegistry_))
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
    Intragranular_bubble_concentration_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intragranular_bubble_concentration", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intragranular_bubble_radius_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intragranular_bubble_radius", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intergranular_bubble_concentration_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_bubble_concentration", 
            sciantixFieldsDict_, 
            addToRegistry_, 
            4.0e+13
        )
    ),
    Intergranular_atoms_per_bubble_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_atoms_per_bubble", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intergranular_vacancies_per_bubble_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_vacancies_per_bubble", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intergranular_bubble_radius_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_bubble_radius", 
            sciantixFieldsDict_, 
            addToRegistry_, 
            1.0e-08
        )
    ),
    Intergranular_bubble_area_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_bubble_area", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intergranular_bubble_volume_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_bubble_volume", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intergranular_fractional_coverage_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_fractional_coverage", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ),
    Intergranular_saturation_fractional_coverage_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_saturation_fractional_coverage", 
            sciantixFieldsDict_, 
            addToRegistry_,
            0.5
        )
    ), 
    Intergranular_fractional_intactness_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::Intergranular_fractional_intactness", 
            sciantixFieldsDict_, 
            addToRegistry_, 
            1.0
        )
    ), 
    timeStepNumber_
    (
        createOrRead<scalar>
        (
            mesh, 
            "SCIANTIX::timeStepNumber", 
            sciantixFieldsDict_, 
            addToRegistry_
        )
    ), 
    grainModes_(0),
    time_(mesh_.time().value()),
    outerIteration_(0),
    nFrequency_
    (       
    #ifdef OPENFOAMFOUNDATION
        fissionGasReleaseDict_.isNull() ? 
    #elif OPENFOAMESI
        fissionGasReleaseDict_.isNullDict() ? 
    #endif             
        1 :
        fissionGasReleaseDict_.lookupOrDefault("nFrequency", 1)
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
    Q_(nullptr),
    Bu_(nullptr),
    T_(nullptr),
    rho_(nullptr),
    sigma_(nullptr),
    inputRead_(false),
    fgrPercent_(0),
    fgrM3_(0),
    fgpM3_(0)
{ 
    // Create grain modes scalar fields
    //  To simplify the costructor these fileds are created as a PtrList.
    //  However, this means that they cannot be references, but local fields.
    //  If addToRegistry is true, the writing of these fields must be done
    //  manually (see function write(os)).

    grainModes_.setSize(20);
    forAll(grainModes_, i)
    {
        grainModes_.set
        (
            i,
            new scalarField
            (
                // addToRegistry must always be false for the grain modes,
                // otherwise two pointers (one in this class, one in globalFields)
                // are created for the same object
                createOrRead<scalar>
                (
                    mesh, 
                    "SCIANTIX::grainModes_" + std::to_string(i),
                    sciantixFieldsDict_, 
                    false
                )
            )
        );
    } 

    // Calculate current fgr value
    scalarField produced((mesh_.V()*Gas_produced_*pos(grainRadius_)));
    scalarField released((mesh_.V()*Gas_released_*pos(grainRadius_)));

    if(refabricatedRod_)
    {   
        if(gSum(Gas_released_t0_) < SMALL)
        {
            Gas_released_t0_ = Gas_released_;
        }

        scalarField released0((mesh_.V()*Gas_released_t0_*pos(grainRadius_)));
        fgr_ = (released - released0)/max(produced - released0, SMALL)*100;
    }
    else
    {
        fgr_ = released/max(produced, SMALL)*100;
    }

    fgr0_ = fgr_;

    // Initialize time, temp, fission rate and sigmaH list to hold 2 values
    Input_history_points = 2;
    Time_input.resize(Input_history_points);
    Temperature_input.resize(Input_history_points);
    Fissionrate_input.resize(Input_history_points);
    Hydrostaticstress_input.resize(Input_history_points);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fgrSCIANTIX::~fgrSCIANTIX()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fgrSCIANTIX::setSCIANTIXVariables
(
        const label cellI
)
{    
    // Reference to Burnup
    const scalarField& BuiOld = Bu_->oldTime().internalField();
    
    //Fuel density
    const scalarField& rho = rho_->internalField(); 
    
    Time_step_number = timeStepNumber_[cellI];

    for(int i=0; i < 20; i++)
    {
        gas_grain_modes[i] = grainModes_[i][cellI];
    }

    Sciantix_variables[0] = grainRadius_[cellI];
    Sciantix_variables[1] = Gas_produced_[cellI];
    Sciantix_variables[2] = Gas_grain_[cellI];
    Sciantix_variables[3] = Gas_grain_solution_[cellI];
    Sciantix_variables[4] = Gas_grain_bubbles_[cellI];
    Sciantix_variables[5] = Gas_boundary_[cellI];
    Sciantix_variables[6] = Gas_released_[cellI]; 

    Sciantix_variables[20] = BuiOld[cellI]/1000;
    Sciantix_variables[21] = Effective_burn_up_[cellI];
    Sciantix_variables[84] = rho[cellI];
    Sciantix_variables[85] = oxygenMetalRatio_[cellI];
    Sciantix_variables[22] = Helium_produced_[cellI];
    Sciantix_variables[86] = Helium_produced_[cellI]; 
    Sciantix_variables[87] = Helium_grain_[cellI];
    Sciantix_variables[88] = Helium_grain_solution_[cellI];
    Sciantix_variables[89] = Helium_grain_bubbles_[cellI];
    Sciantix_variables[90] = Helium_boundary_[cellI];
    Sciantix_variables[91] = Helium_released_[cellI];

    Sciantix_variables[7] = Intragranular_bubble_concentration_[cellI];
    Sciantix_variables[8] = Intragranular_bubble_radius_[cellI];
    Sciantix_variables[9] = tr(intragranularGasSwelling_[cellI])/3;
    Sciantix_variables[10] = Intergranular_bubble_concentration_[cellI];
    Sciantix_variables[11] = Intergranular_atoms_per_bubble_[cellI];
    Sciantix_variables[12] = Intergranular_vacancies_per_bubble_[cellI];
    Sciantix_variables[13] = Intergranular_bubble_radius_[cellI];
    Sciantix_variables[14] = Intergranular_bubble_area_[cellI];
    Sciantix_variables[15] = Intergranular_bubble_volume_[cellI];
    Sciantix_variables[16] = Intergranular_fractional_coverage_[cellI];
    Sciantix_variables[17] = Intergranular_saturation_fractional_coverage_[cellI];
    Sciantix_variables[18] = tr(intergranularGasSwelling_[cellI])/3;
    Sciantix_variables[19] = Intergranular_fractional_intactness_[cellI]; 
}


void Foam::fgrSCIANTIX::updateSCIANTIXVariables
(
        const label cellI
)
{    
    timeStepNumber_[cellI] = Time_step_number;

    for(int i=0; i < 20; i++)
    {
        grainModes_[i][cellI] = gas_grain_modes[i];
    }

    grainRadius_[cellI] = Sciantix_variables[0];
    Gas_produced_[cellI] = Sciantix_variables[1];
    Gas_grain_[cellI] = Sciantix_variables[2];
    Gas_grain_solution_[cellI] = Sciantix_variables[3];
    Gas_grain_bubbles_[cellI] = Sciantix_variables[4];
    Gas_boundary_[cellI] = Sciantix_variables[5];
    Gas_released_[cellI] = Sciantix_variables[6];

    Effective_burn_up_[cellI] = Sciantix_variables[21];
    oxygenMetalRatio_[cellI] = Sciantix_variables[85];

    Helium_produced_[cellI] = Sciantix_variables[86];
    Helium_grain_[cellI] = Sciantix_variables[87];
    Helium_grain_solution_[cellI] = Sciantix_variables[88];
    Helium_grain_bubbles_[cellI] = Sciantix_variables[89];
    Helium_boundary_[cellI] = Sciantix_variables[90];
    Helium_released_[cellI] = Sciantix_variables[91];

    Intragranular_bubble_concentration_[cellI] = Sciantix_variables[7];
    Intragranular_bubble_radius_[cellI] = Sciantix_variables[8];
    intragranularGasSwelling_[cellI] = Sciantix_variables[9]*I;
    Intergranular_bubble_concentration_[cellI] = Sciantix_variables[10];
    Intergranular_atoms_per_bubble_[cellI] = Sciantix_variables[11];
    Intergranular_vacancies_per_bubble_[cellI] = Sciantix_variables[12];
    Intergranular_bubble_radius_[cellI] = Sciantix_variables[13];
    Intergranular_bubble_area_[cellI] = Sciantix_variables[14];
    Intergranular_bubble_volume_[cellI] = Sciantix_variables[15];
    Intergranular_fractional_coverage_[cellI] = Sciantix_variables[16];
    Intergranular_saturation_fractional_coverage_[cellI] = Sciantix_variables[17];
    intergranularGasSwelling_[cellI] = Sciantix_variables[18]*I;
    Intergranular_fractional_intactness_[cellI] = Sciantix_variables[19];

}


void Foam::fgrSCIANTIX::correct()
{  
    if
    (
        Q_ == nullptr or Bu_ == nullptr or T_ == nullptr or rho_ == nullptr
        or sigma_ == nullptr
    )
    {
        Q_ = &mesh_.lookupObject<volScalarField>("Q");
        Bu_ = &mesh_.lookupObject<volScalarField>("Bu");
        T_ = &mesh_.lookupObject<volScalarField>("T");
        rho_ = &mesh_.lookupObject<volScalarField>("rho");
        sigma_ = &mesh_.lookupObject<volSymmTensorField>("sigma");
    }

    intragranularGasSwelling_.storeOldTimes();
    intergranularGasSwelling_.storeOldTimes();

    if(mesh_.time().value() > time_)
    {
        Time_input[0] = time_/3600;
        Time_input[1] = mesh_.time().value()/3600;
        time_ = mesh_.time().value();

        Time_end_h = Time_input[Input_history_points-1];
        Time_end_s = Time_end_h * 3600.0;

        if(!inputRead_)
        {
            if( inputSCIANTIX_.size() > 0 )
            {
                // Call OFFBEAT's input reading function (from dict)
                readSCIANTIXinput();
            }
            else
            {
                // Call SCIANTIX's input reading function (from file)
                InputReading();
            }
            
            SetOptions( Sciantix_options,
                Sciantix_scaling_factors );

            inputRead_ = true;
        }

        outerIteration_ = 0;
    }

    if(outerIteration_ % nFrequency_ == 0)
    {
        const scalarField& Gas_releasedOld(Gas_released_);
        const scalarField& Helium_releasedOld(Helium_released_);

        const scalarField& Q = Q_->internalField();
        const scalarField& QOld = Q_->oldTime().internalField();

        // Hydrostatic stress
        const symmTensorField& sigma = sigma_->internalField();
        const symmTensorField& sigmaOld = sigma_->oldTime().internalField();

        const scalarField sigmaH = tr(sigma)/3/1e6;  
        const scalarField sigmaHOld = tr(sigmaOld)/3/1e6;  

        const scalarField& Ti = T_->internalField();
        const scalarField& TiOld = T_->oldTime().internalField();

        // Fuel cell volumes [m^3]
        const scalarField& V = mesh_.V();

        // Put the number of helium and gas moles released to zero
        scalar molFGRprev_(molFGR_);
        molHe_ *= 0;
        molFGR_ *= 0;

        const labelListList& matAddrList(mat_.matAddrList());

        forAll(mat_.materialsList(), i)
        {
            if(isA<fuelMaterial>(mat_.materialsList()[i]))
            {
                const labelList& addr(matAddrList[i]);

                forAll(addr, j)
                {
                    const label cellI = addr[j];

                    setSCIANTIXVariables(cellI);

                    Temperature_input[0] = TiOld[cellI];
                    Temperature_input[1] = Ti[cellI];

                    Fissionrate_input[0] = max(QOld[cellI]/312.0e-13, SMALL);
                    Fissionrate_input[1] = Q[cellI]/312.0e-13;

                    Hydrostaticstress_input[0] = sigmaHOld[cellI];
                    Hydrostaticstress_input[1] = sigmaH[cellI];
                
                    dTime_s = (Time_input[1] - Time_input[0])*3600;
                    dTime_h = Time_input[1] - Time_input[0];

                    RunSCIANTIX( );    

                    intragranularGasSwelling_[cellI] = 
                    (
                        relax_*Sciantix_variables[9]*I 
                        + (1- relax_)*intragranularGasSwelling_[cellI]
                    );

                    intergranularGasSwelling_[cellI] = 
                    (
                        relax_*Sciantix_variables[18]*I 
                        + (1- relax_)*intergranularGasSwelling_[cellI]
                    );

                    //***** Helium released *****//
                    molHe_ += 
                    (
                        Sciantix_variables[91] 
                        - Helium_releasedOld[cellI]
                    )*V[cellI]/6.02e23;  

                    //**** Total gas released *****//
                    molFGR_ += 
                    (
                        Sciantix_variables[6] 
                        - Gas_releasedOld[cellI]
                    )*V[cellI]/6.02e23;   

                }
            }
        }

        reduce(molFGR_, sumOp<scalar>());

        molFGR_ = (1- relax_)*molFGRprev_ + relax_*molFGR_;
    }

    outerIteration_++;
}

void Foam::fgrSCIANTIX::updateVariables()
{    
    if(mesh_.time().value() > time_)
    {
        Time_input[0] = time_/3600;
        Time_input[1] = mesh_.time().value()/3600;
        time_ = mesh_.time().value();

        Time_end_h = Time_input[Input_history_points-1];
        Time_end_s = Time_end_h * 3600.0;

        if(!inputRead_)
        {
            if( inputSCIANTIX_.size() > 0  )
            {
                // Call OFFBEAT's input reading function (from dict)
                readSCIANTIXinput();
            }
            else
            {
                // Call SCIANTIX's input reading function (from file)
                InputReading();
            }
            
            SetOptions( Sciantix_options,
                Sciantix_scaling_factors );

            inputRead_ = true;
        }
    }

    const scalarField& Q = Q_->internalField();
    const scalarField& QOld = Q_->oldTime().internalField();

    //Hydrostatic stress
    const symmTensorField& sigma = sigma_->internalField();
    const symmTensorField& sigmaOld = sigma_->oldTime().internalField();

    const scalarField sigmaH = tr(sigma)/3/1e6;  
    const scalarField sigmaHOld = tr(sigmaOld)/3/1e6;  

    Input_history_points = 2;

    const scalarField& Ti = T_->internalField();
    const scalarField& TiOld = T_->oldTime().internalField();

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
       if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];
                setSCIANTIXVariables(cellI);

                Temperature_input[0] = TiOld[cellI];
                Temperature_input[1] = Ti[cellI];

                Fissionrate_input[0] = max(QOld[cellI]/312.0e-13, SMALL);
                Fissionrate_input[1] = Q[cellI]/312.0e-13;

                Hydrostaticstress_input[0] = sigmaHOld[cellI];
                Hydrostaticstress_input[1] = sigmaH[cellI];
                
                dTime_s = (Time_input[1] - Time_input[0])*3600;
                dTime_h = Time_input[1] - Time_input[0];

                RunSCIANTIX( );

                updateSCIANTIXVariables(cellI);
            }
        }
    }

    scalarField produced(mesh_.V()*pos(grainRadius_)*Gas_produced_);
    scalarField released(mesh_.V()*pos(grainRadius_)*Gas_released_);
    scalarField released0((mesh_.V()*Gas_released_t0_*pos(grainRadius_)));

    // Update fgr and fgr0 (for time stepping)
    fgr0_ = fgr_;
    
    if(refabricatedRod_)
    {
        fgr_ = (released - released0)/max(produced - released0, SMALL)*100;
    }
    else
    {
        fgr_ = released/max(produced,SMALL)*100;
    }

    // scalar HeProduced(gSum(mesh_.V()*Helium_produced_*pos(grainRadius_)));
    // scalar HeReleased(gSum(mesh_.V()*Helium_released_*pos(grainRadius_)));

    // For debugging
    Info << "fission gas (Xe/Kr) atoms produced " << gSum(produced) << endl;
    Info << "fission gas (Xe/Kr) atoms released " << gSum(released) << endl;

    const globalOptions& globalOpt
    (
        mesh_.lookupObject<globalOptions>("globalOptions")
    );

    // --- Update quantities for post processing
    
    // Update fgrM3 (for functionObject)
    const scalar fgrAtoms = 
    gSum(mesh_.V()*Gas_released_*pos(grainRadius_))/globalOpt.angularFraction();

    // Use ideal gas law FGR (m3) @ T= 293 K and P=101325 Pa
    scalar nMolesReleased = fgrAtoms/6.022e23; 
    fgrM3_  = nMolesReleased*8.314*293/101325; 

    // Update fgrPercent (for functionObject)
    fgrPercent_ = gSum(released)/max(gSum(produced),SMALL)*100;

    // Update fgpM3 (for functionObject)
    const scalar fgpAtoms = 
    gSum(mesh_.V()*Gas_produced_*pos(grainRadius_)/globalOpt.angularFraction());

    // Use ideal gas law fgp (m3) @ T= 293 K and P=101325 Pa
    scalar nMolesProduced = fgpAtoms/6.022e23; 
    fgpM3_  = nMolesProduced*8.314*293/101325; 
    
    // --- End update quantities for post processing
    
    if(refabricatedRod_)
    {
        Info << "fgr% " << fgrPercent_ << endl;
        
        Info << "fgr% refabricated " << 
            gSum(released - released0)/max(
                gSum(produced - released0),SMALL)*100 << nl << endl;
    }
    else
    {
        Info << "fgr% " << fgrPercent_ << endl;
    }


    // Info << "He atoms produced " << HeProduced << endl;
    // Info << "He atoms released " << HeReleased << endl;
    // Info << "He% " << HeReleased/max(HeProduced,SMALL)*100 << nl << endl;

    grainRadius_.correctBoundaryConditions();
    oxygenMetalRatio_.correctBoundaryConditions();
    intergranularGasSwelling_.correctBoundaryConditions();
    intragranularGasSwelling_.correctBoundaryConditions();
}


const Foam::wordList Foam::fgrSCIANTIX::gasComponents() const
{
        wordList gasComp(3);
        gasComp[0] = "Xe";
        gasComp[1] = "Kr";
        gasComp[2] = "He";

        return gasComp;
}


const Foam::scalarList Foam::fgrSCIANTIX::gasMols() const
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


bool Foam::fgrSCIANTIX::writeData(Ostream& os) const
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

    if(refabricatedRod_)
    {        
        scalar deltaProduced(gSum(mesh_.V()*(
            Gas_produced_ -  Gas_released_t0_)*pos(grainRadius_)));
        
        scalar deltaReleased(gSum(mesh_.V()*(
            Gas_released_ - Gas_released_t0_)*pos(grainRadius_)));

        Info << "fgr% refabricated " << 
            deltaReleased/max(deltaProduced,SMALL)*100 << nl << endl;
    }

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
        
        if(refabricatedRod_)
        {
            os.writeKeyword("SCIANTIX::Gas_released_t0") 
            << Gas_released_t0_ << token::END_STATEMENT << nl << nl;
        }

        os.writeKeyword("SCIANTIX::Effective_burn_up") 
        << Effective_burn_up_ << token::END_STATEMENT << nl << nl;
        // os.writeKeyword("Oxygen_to_metal_ratio") 
        // << oxygenMetalRatio_ << token::END_STATEMENT << nl << nl;
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
        os.writeKeyword("SCIANTIX::Intragranular_bubble_concentration") 
        << Intragranular_bubble_concentration_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intragranular_bubble_radius") 
        << Intragranular_bubble_radius_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_bubble_concentration") 
        << Intergranular_bubble_concentration_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_atoms_per_bubble") 
        << Intergranular_atoms_per_bubble_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_vacancies_per_bubble") 
        << Intergranular_vacancies_per_bubble_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_bubble_radius") 
        << Intergranular_bubble_radius_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_bubble_area") 
        << Intergranular_bubble_area_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_bubble_volume") 
        << Intergranular_bubble_volume_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_fractional_coverage") 
        << Intergranular_fractional_coverage_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_saturation_fractional_coverage") 
        << Intergranular_saturation_fractional_coverage_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::Intergranular_fractional_intactness") 
        << Intergranular_fractional_intactness_ << token::END_STATEMENT << nl << nl;
        os.writeKeyword("SCIANTIX::timeStepNumber") 
        << timeStepNumber_ << token::END_STATEMENT << nl << nl;

        forAll(grainModes_, i)
        {
            word fieldName("SCIANTIX::grainModes_" + std::to_string(i));
            os.writeKeyword(fieldName)
            << grainModes_[i] << token::END_STATEMENT << nl << nl;
        } 
    }
    else
    {
        //  It is necessary to take care separately of the grainModes fields
        //  as the fact that they are in a list implies that they cannot be 
        //  references but local fields. Thus, the internal field of the 
        //  corresponding volumetric field (if addToRegistry is used)
        //  is not automatically updated.

        //Create grainMode field
        forAll(grainModes_, i)
        {
            tmp<volScalarField> tgrainModesField
            (
                new volScalarField
                (
                    IOobject
                    (
                        "SCIANTIX::grainModes_" + std::to_string(i),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "SCIANTIX::grainModes_" + std::to_string(i),
                        dimless,
                        0.0
                    ),
                    "zeroGradient"
                )

            );
            scalarField& grainModesI(tgrainModesField->ref());
            grainModesI = grainModes_[i];            
            tgrainModesField->correctBoundaryConditions();
            tgrainModesField->write();
        } 
    }
    
    return os.good();
}
// ************************************************************************* //
