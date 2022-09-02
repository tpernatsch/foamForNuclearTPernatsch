model test
  extends OpenModelica;
  Modelica.Blocks.Sources.Ramp ramp(duration = 15, height = 1, offset = 0, startTime = 5)  annotation(
    Placement(visible = true, transformation(origin = {-10, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput momentumModelica annotation(
    Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput power annotation(
    Placement(visible = true, transformation(origin = {-28, -42}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-28, -42}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator integrator annotation(
    Placement(visible = true, transformation(origin = {38, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput powerIntegral annotation(
    Placement(visible = true, transformation(origin = {82, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(ramp.y, momentumModelica) annotation(
    Line(points = {{2, -2}, {50, -2}, {50, 0}}, color = {0, 0, 127}));
  connect(power, integrator.u) annotation(
    Line(points = {{-28, -42}, {3, -42}, {3, -46}, {26, -46}}, color = {0, 0, 127}));
  connect(integrator.y, powerIntegral) annotation(
    Line(points = {{50, -46}, {82, -46}, {82, -50}}, color = {0, 0, 127}));

annotation(
    uses(Modelica(version = "3.2.3")));
end test;