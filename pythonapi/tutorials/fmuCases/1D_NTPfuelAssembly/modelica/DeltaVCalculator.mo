model DeltaVCalculator
  Modelica.Blocks.Continuous.Integrator integrator annotation(
    Placement(transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput thrust annotation(
    Placement(transformation(origin = {-120, 6}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput deltaV annotation(
    Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  parameter Modelica.SIunits.Mass dryMass(min = 1) = 20000 "Dry mass" annotation(
    Placement(transformation(origin = {-90, -10}, extent = {{-10, -10}, {10, 10}})));
  parameter Modelica.SIunits.Mass propellantMass(min = 1) = 20000 "Propellant mass" annotation(
    Placement(transformation(origin = {-90, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Division division annotation(
    Placement(transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput massFlow annotation(
    Placement(transformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Continuous.Integrator propellantIntegrator annotation(
    Placement(transformation(origin = {-70, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.GreaterEqualThreshold greaterEqualThreshold(threshold = propellantMass) annotation(
    Placement(transformation(origin = {-10, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interaction.Show.BooleanValue endPropellant annotation(
    Placement(transformation(origin = {30, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add totalMass annotation(
    Placement(transformation(origin = {10, -26}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add totalMassPropellant(k2 = -1) annotation(
    Placement(transformation(origin = {-30, -50}, extent = {{-10, -10}, {10, 10}})));
equation
  totalMass.u1 = dryMass;
  totalMassPropellant.u1 = propellantMass;
  totalMassPropellant.u2 = min(propellantIntegrator.y, propellantMass);
  connect(thrust, division.u1) annotation(
    Line(points = {{-120, 6}, {38, 6}}, color = {0, 0, 127}));
  connect(division.y, integrator.u) annotation(
    Line(points = {{61, 0}, {68, 0}}, color = {0, 0, 127}));
  connect(integrator.y, deltaV) annotation(
    Line(points = {{91, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(massFlow, propellantIntegrator.u) annotation(
    Line(points = {{-120, -80}, {-82, -80}}, color = {0, 0, 127}));
  connect(propellantIntegrator.y, greaterEqualThreshold.u) annotation(
    Line(points = {{-58, -80}, {-22, -80}}, color = {0, 0, 127}));
  connect(greaterEqualThreshold.y, endPropellant.activePort) annotation(
    Line(points = {{2, -80}, {18, -80}}, color = {255, 0, 255}));
  connect(totalMass.y, division.u2) annotation(
    Line(points = {{21, -26}, {27, -26}, {27, -6}, {37, -6}}, color = {0, 0, 127}));
  connect(totalMassPropellant.y, totalMass.u2) annotation(
    Line(points = {{-18, -50}, {-10, -50}, {-10, -32}, {-2, -32}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end DeltaVCalculator;