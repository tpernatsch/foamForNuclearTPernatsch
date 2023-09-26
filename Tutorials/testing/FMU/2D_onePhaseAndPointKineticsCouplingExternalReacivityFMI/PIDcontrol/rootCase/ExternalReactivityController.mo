model ExternalReactivityController
  Modelica.Blocks.Interfaces.RealOutput externalReactivity(start = 0) annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput power(start = 10e6) annotation(
    Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 1e6, offset = 10e6, startTime = 110) annotation(
    Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Ziegler–Nichols tuning method (no overshoot)
  parameter Real Ku = 2e-11;
  parameter Real Tu = 5;
  parameter Real Kp = 0.2 * Ku;
  parameter Real Ti = 0.5 * Tu;
  parameter Real Td = 0.33 * Tu;
  Modelica.Blocks.Continuous.LimPID pid(k = Kp, Ti = Ti, Td = Td, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.InitPID.NoInit, limitsAtInit = true, yMax = 50e-5) annotation(
    Placement(visible = true, transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  /*1e-12*/
  /*0.1*/
  Modelica.Blocks.Logical.Switch switch1 annotation(
    Placement(visible = true, transformation(origin = {-20, 8}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = 0.5) annotation(
    Placement(visible = true, transformation(origin = {-60, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step step(height = 1, offset = 0, startTime = 20) annotation(
    Placement(visible = true, transformation(origin = {-90, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(ramp.y, pid.u_s) annotation(
    Line(points = {{-59, 50}, {-12, 50}}, color = {0, 0, 127}));
  connect(pid.y, externalReactivity) annotation(
    Line(points = {{12, 50}, {60, 50}, {60, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(switch1.y, pid.u_m) annotation(
    Line(points = {{-9, 8}, {0, 8}, {0, 38}}, color = {0, 0, 127}));
  connect(power, switch1.u1) annotation(
    Line(points = {{-120, 0}, {-32, 0}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, switch1.u2) annotation(
    Line(points = {{-49, 20}, {-39.5, 20}, {-39.5, 8}, {-32, 8}}, color = {255, 0, 255}));
  connect(step.y, greaterThreshold.u) annotation(
    Line(points = {{-79, 20}, {-73, 20}}, color = {0, 0, 127}));
  connect(ramp.y, switch1.u3) annotation(
    Line(points = {{-58, 50}, {-36, 50}, {-36, 16}, {-32, 16}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end ExternalReactivityController;