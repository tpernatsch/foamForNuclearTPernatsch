model Controller
  Modelica.Blocks.Interfaces.RealOutput pumpMomentum annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp massFlowCmd(duration = 10, offset = 590, startTime = 100, height = -300) annotation(
    Placement(transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID pid(k = 0.0001, Ti = 0.1, controllerType = Modelica.Blocks.Types.SimpleController.PI, limitsAtInit = true, yMin = 0) annotation(
    Placement(transformation(extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput massFlow annotation(
    Placement(transformation(origin = {-120, -30}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = 2)  annotation(
    Placement(transformation(origin = {-30, -30}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(massFlowCmd.y, pid.u_s) annotation(
    Line(points = {{-59, 0}, {-12, 0}}, color = {0, 0, 127}, thickness = 0.5));
  connect(pid.y, pumpMomentum) annotation(
    Line(points = {{11, 0}, {110, 0}}, color = {0, 0, 127}, thickness = 0.5));
  connect(firstOrder.y, pid.u_m) annotation(
    Line(points = {{-18, -30}, {0, -30}, {0, -12}}, color = {0, 0, 127}));
  connect(massFlow, firstOrder.u) annotation(
    Line(points = {{-120, -30}, {-42, -30}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end Controller;