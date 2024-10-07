///////////////////////////////////////////
//                                       //
//           S C I A N T I X             //
//           ---------------             //
//                                       //
//  Version: 1.4                         //
//  Year   : 2019                        //
//  Authors: D. Pizzocri and T. Barani   //
//           L. Cognini                  //
///////////////////////////////////////////

#include "SetOptions.h"
#include <iostream>

void SetOptions(int Sciantix_options[], double Sciantix_scaling_factors[])
{
  iverification = Sciantix_options[0];
  igrain_growth = Sciantix_options[1];
  iinert_gas_behavior = Sciantix_options[2];
  igas_diffusion_coefficient = Sciantix_options[3];
  iintra_bubble_evolution = Sciantix_options[4];
  ibubble_radius = Sciantix_options[5];
  iresolution_rate = Sciantix_options[6];
  itrapping_rate = Sciantix_options[7];
  inucleation_rate = Sciantix_options[8];
  isolver = Sciantix_options[9];
  ioutput = Sciantix_options[10];
  igrain_boundary_vacancy_diffusion_coefficient = Sciantix_options[11];
  igrain_boundary_behaviour = Sciantix_options[12];
  igrain_boundary_micro_cracking = Sciantix_options[13];
  igrain_recrystallization = Sciantix_options[14];
  ifuel_reactor_type = Sciantix_options[15];
  igas_effective_coefficient = Sciantix_options[16];
  igas_sweeping = Sciantix_options[17];
  imicro_cracking_span = Sciantix_options[18];

  sf_resolution_rate = Sciantix_scaling_factors[0];
  sf_trapping_rate = Sciantix_scaling_factors[1];
  sf_nucleation_rate = Sciantix_scaling_factors[2];
  sf_diffusion_rate = Sciantix_scaling_factors[3];
}
