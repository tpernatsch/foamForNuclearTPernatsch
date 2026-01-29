///////////////////////////////////////////
//                                       //
//           S C I A N T I X             //
//           ---------------             //
//                                       //
//  Version: 1.4                         //
//  Year   : 2019                        //
//  Authors: D. Pizzocri and T. Barani   //
//                                       //
///////////////////////////////////////////

/// GasDiffusion
/// This function calls the solver SpectralDiffusion
/// to evaluate the fission gas concentration
/// inside the single grain after a fixed time step.
/// Typical numerical values for UO2 fuel are taken from:
/// [1] Olander, Wongsawaeng, Journal of Nuclear Materials, 354 (2006), 94-109

#include "GasDiffusion.h"
#include <iostream>

void GasDiffusion( )
{
  double diffusion_coefficient = GasDiffusionCoefficient(Temperature[0], Fissionrate[0]);
  double resolution_rate = ResolutionRate(Intragranular_bubble_radius[0], Fissionrate[0]);

  double trapping_rate = TrappingRate(diffusion_coefficient, Intragranular_bubble_radius[0], Intragranular_bubble_concentration[0]);
  double equilibrium_fraction(resolution_rate / (resolution_rate + trapping_rate));

  switch(igas_effective_coefficient)
  {
    case 0 :
      // Correlation provides effective coefficient
      diffusion_coefficient *= 1;
      break;

    case 1 :
      // Correlation provides single atom coefficient
      diffusion_coefficient *= equilibrium_fraction;      
      break;

    default :
      ErrorMessages::Switch("EffectiveDiffusionCoefficient", "igas_effective_coefficient", igas_effective_coefficient);
      break;
  }

  const unsigned short int N(20);

  const double fission_yield = Fission_yield_Xe + Fission_yield_Kr;
  double source_term = fission_yield * Fissionrate[1]; // (at/m3-s)
  double initial_condition_term[4] = {gas_grain_modes[0], gas_grain_modes[1], gas_grain_modes[2], gas_grain_modes[3]};  

  switch(isolver)
  {
    case 0 :
      Gas_grain[1] = Solver::FORMAS(initial_condition_term, diffusion_coefficient, Grain_radius[1], source_term, dTime_s);          
      
      gas_grain_modes[0] = initial_condition_term[0];
      gas_grain_modes[1] = initial_condition_term[1];
      gas_grain_modes[2] = initial_condition_term[2];
      gas_grain_modes[3] = initial_condition_term[3];
      break;

    case 1 :
      Gas_grain[1] = Solver::SpectralDiffusion(gas_grain_modes, N, diffusion_coefficient, Grain_radius[1], source_term, dTime_s);
      break;

    default :
      ErrorMessages::Switch("GasDiffusion", "isolver", isolver);
      break;
  }

  Gas_grain_solution[1] = Gas_grain[1] * equilibrium_fraction;
  Gas_grain_bubbles[1] = Gas_grain[1] * (1.0 - equilibrium_fraction);

  Gas_boundary[1] = Gas_produced[1] - Gas_grain[1] - Gas_released[1];

  //- Boundary Sweeping
  switch(igas_sweeping)
  {
    case 0 :
      // do nothing
      break;

    case 1 :
    {
      double oldVolume = pow(Grain_radius[0], 3.0);
      double newVolume = pow(Grain_radius[1], 3.0);
      double VolumDiff = (newVolume- oldVolume)/newVolume;
      Gas_grain[1] -= VolumDiff * Gas_grain[1];
      Gas_boundary[1] += VolumDiff * Gas_grain[1];
      for(int i=0; i<20; i++)
      {
        gas_grain_modes[i] *= (1 - VolumDiff);
      }  

      break;
    }

    default :
      ErrorMessages::Switch("GasDiffusion", "igas_sweeping", igas_sweeping);
      break;
  }
}
