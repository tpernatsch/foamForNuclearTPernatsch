model lumpedSystem
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = 0.1, initType = Modelica.Blocks.Types.Init.InitialOutput, k = 0.05, y(start = 0))  annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput u annotation(
    Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-68, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput y annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = 50)  annotation(
    Placement(visible = true, transformation(origin = {0, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step step(height = -0.01, startTime = 2.5) annotation(
    Placement(visible = true, transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add3 add3 annotation(
    Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(firstOrder.u, u) annotation(
    Line(points = {{-12, 0}, {-120, 0}}, color = {0, 0, 127}));
  connect(add3.y, y) annotation(
    Line(points = {{72, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(step.y, add3.u1) annotation(
    Line(points = {{12, 50}, {40, 50}, {40, 8}, {48, 8}}, color = {0, 0, 127}));
  connect(const.y, add3.u3) annotation(
    Line(points = {{12, -50}, {40, -50}, {40, -8}, {48, -8}}, color = {0, 0, 127}));
  connect(firstOrder.y, add3.u2) annotation(
    Line(points = {{12, 0}, {48, 0}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end lumpedSystem;