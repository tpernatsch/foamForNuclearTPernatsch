model momentumSourceTest
  // extends OpenModelica;
  Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 1, offset = 1, startTime = 20) annotation(
    Placement(transformation(origin = {-50, -20}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput modMomentumSourceOut1 annotation(
    Placement(transformation(origin = {110, -20}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput modTin annotation(
    Placement(visible = true, transformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 50, offset = 350, startTime = 20) annotation(
    Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput modTemperatureOut annotation(
    Placement(transformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput modHout1 annotation(
    Placement(transformation(origin = {110, 80}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp2(duration = 10, height = -5000, offset = 10000, startTime = 20) annotation(
    Placement(visible = true, transformation(origin = {-50, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput modPowerIn annotation(
    Placement(visible = true, transformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator1 annotation(
    Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput hxPowerIntegral annotation(
    Placement(transformation(origin = {110, -80}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput modHout2 annotation(
    Placement(transformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput modMomentumSourceOut2 annotation(
    Placement(transformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp11(duration = 10, height = 61.6e6, offset = 50e6, startTime = 20) annotation(
    Placement(transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput modPowerOut annotation(
    Placement(transformation(origin = {110, 10}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(ramp.y, modMomentumSourceOut1) annotation(
    Line(points = {{-39, -20}, {110, -20}}, color = {0, 0, 127}));
  connect(ramp1.y, modTemperatureOut) annotation(
    Line(points = {{-39, 40}, {110, 40}}, color = {0, 0, 127}));
  connect(ramp2.y, modHout1) annotation(
    Line(points = {{-39, 80}, {110, 80}}, color = {0, 0, 127}));
  connect(modPowerIn, integrator1.u) annotation(
    Line(points = {{-120, -80}, {-12, -80}}, color = {0, 0, 127}));
  connect(integrator1.y, hxPowerIntegral) annotation(
    Line(points = {{12, -80}, {110, -80}}, color = {0, 0, 127}));
  connect(ramp2.y, modHout2) annotation(
    Line(points = {{-38, 80}, {80, 80}, {80, 60}, {110, 60}}, color = {0, 0, 127}));
  connect(ramp11.y, modPowerOut) annotation(
    Line(points = {{-39, 10}, {110, 10}}, color = {0, 0, 127}));
  connect(ramp.y, modMomentumSourceOut2) annotation(
    Line(points = {{-38, -20}, {80, -20}, {80, -40}, {110, -40}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end momentumSourceTest;