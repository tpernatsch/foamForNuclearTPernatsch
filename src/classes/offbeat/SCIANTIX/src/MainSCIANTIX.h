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

#include "MainVariables.h"
#include "Sciantix.h"
#include "SolverVerification.h"
#include "Solver.h"
#include "InputInterpolation.h"
#include "InputReading.h"
#include "SetOptions.h"
#include "Initialization.h"
#include "TimeStepCalculation.h"
#include "OutputWriting.h"

void RunSCIANTIX( );
