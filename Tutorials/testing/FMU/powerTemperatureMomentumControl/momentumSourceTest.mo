model momentumSourceTest
  // extends OpenModelica;
  Modelica.Blocks.Sources.Ramp ramp(duration = 15, height = 1, offset = 0, startTime = 5)  annotation(
    Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput momentumModelica annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput power annotation(
    Placement(visible = true, transformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator annotation(
    Placement(visible = true, transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput powerIntegral annotation(
    Placement(visible = true, transformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp1(duration = 15, height = 100, offset = 850, startTime = 5) annotation(
    Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput temperatureModelica annotation(
    Placement(visible = true, transformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput HModelica annotation(
    Placement(visible = true, transformation(origin = {110, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp2(duration = 5, height = -9990, offset = 10000, startTime = 15) annotation(
    Placement(visible = true, transformation(origin = {-50, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput hxPower annotation(
    Placement(visible = true, transformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator1 annotation(
    Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput hxPowerIntegral annotation(
    Placement(visible = true, transformation(origin = {110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(ramp.y, momentumModelica) annotation(
    Line(points = {{-39, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(power, integrator.u) annotation(
    Line(points = {{-120, -40}, {-12, -40}}, color = {0, 0, 127}));
  connect(ramp1.y, temperatureModelica) annotation(
    Line(points = {{-39, 40}, {110, 40}}, color = {0, 0, 127}));
  connect(ramp2.y, HModelica) annotation(
    Line(points = {{-39, 80}, {110, 80}}, color = {0, 0, 127}));
  connect(hxPower, integrator1.u) annotation(
    Line(points = {{-120, -80}, {-12, -80}}, color = {0, 0, 127}));
  connect(integrator.y, powerIntegral) annotation(
    Line(points = {{12, -40}, {110, -40}}, color = {0, 0, 127}));
  connect(integrator1.y, hxPowerIntegral) annotation(
    Line(points = {{12, -80}, {110, -80}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end momentumSourceTest;