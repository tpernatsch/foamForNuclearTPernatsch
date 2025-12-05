model ThrustCalculator
  Real M_H2 = 2.014e-3;
  Real gamma_H2 = 1.4;
  Real Rconst = 8.314;
  Real gConst = 9.8;
  Modelica.Blocks.Interfaces.RealInput Tnozzle annotation(
    Placement(transformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput exhaustVelocity annotation(
    Placement(transformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput thrust annotation(
    Placement(transformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput massFlow annotation(
    Placement(transformation(origin = {-120, 60}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 60}, extent = {{-20, -20}, {20, 20}})));
equation
  exhaustVelocity = sqrt(2*gamma_H2/(gamma_H2-1) * Rconst*Tnozzle/M_H2);
  thrust = massFlow * exhaustVelocity;

annotation(
    uses(Modelica(version = "3.2.3")));
end ThrustCalculator;