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

/// SCIANTIX is a 0D code developed at Politecnico di Milano.
/// The objective of SCIANTIX is to represent the behaviour of a single grain of nuclear fuel.
/// The modelling of inert gas behaviour is the main aspect considered.
/// Engineering models are used, allowing for future integration in industrial fuel performance codes.
/// Nevertheless, physically-based model are preferred to empirical models.
/// This facilitates the incorporation of information from lower length scale calculations.

#include "MainSCIANTIX.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>

void RunSCIANTIX( )
{
  //NOTE: these functions take too much time if run for every cell at every iteration.
  //NOTE: they were moved in fgrSCIANTIX.C
  // InputReading( );
  // SetOptions( Sciantix_options,
  //             Sciantix_scaling_factors );
  Initialization( );

  if (Sciantix_options[0]) SolverVerification( );

  //NOTE: these lines take too much computation time. Do we need them in the loop?
  // Output_file.open("output.txt", std::ios::out);
  // Execution_file.open("execution.txt", std::ios::out);

  timer = clock( );

  Time_h = Time_input[0];
  // int times(0);
  while (Time_h <= Time_end_h)
  {
    // std::cout << times << " time, " << Sciantix_history[6] << "\n";
    // times++;
    
    dTime_h = TimeStepCalculation( );
    Sciantix_history[6] = dTime_h * s_h;

    // if (Sciantix_options[10]) OutputWriting( );

    if (Time_h < Time_end_h)
    {
      Time_step_number++;
      Time_h += dTime_h;
      Time_s += Sciantix_history[6];
    }
    else break;

    // Operations to set up the history
    Sciantix_history[0] = Sciantix_history[1];
    Sciantix_history[1] = InputInterpolation(Time_h, Time_input, Temperature_input, Input_history_points);
    Sciantix_history[2] = Sciantix_history[3];
    Sciantix_history[3] = InputInterpolation(Time_h, Time_input, Fissionrate_input, Input_history_points);
    Sciantix_history[4] = Sciantix_history[5];
    Sciantix_history[5] = InputInterpolation(Time_h, Time_input, Hydrostaticstress_input, Input_history_points);
    
    Sciantix(Sciantix_options, Sciantix_history, Sciantix_variables, Sciantix_scaling_factors);

  }

  timer = clock( ) - timer;

  //NOTE: these lines take too much computation time. Do we need them in the loop?
  // Execution_file << std::setprecision(12) << std::scientific << (double)timer/CLOCKS_PER_SEC << "\t" << CLOCKS_PER_SEC << "\t" << (double)timer << "\t" << Time_step_number << std::endl;
  // Output_file.close( );
  // Execution_file.close( );
}
