model ExternalReactivityControllerSimulator
  ExternalReactivityController externalReactivityController annotation(
    Placement(visible = true, transformation(origin = {0, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = 0.5) annotation(
    Placement(visible = true, transformation(origin = {-28, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch1 annotation(
    Placement(visible = true, transformation(origin = {40, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step step(height = -1, offset = 1, startTime = 100) annotation(
    Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp(duration = 2,height = -1e6, offset = 11000000.000000002) annotation(
    Placement(visible = true, transformation(origin = {-70, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  PointKinetics pointKinetics annotation(
    Placement(visible = true, transformation(origin = {0, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  parameter Real Ku = 2e-11;
  parameter Real Tu = 5;
  parameter Real Kp = 0.2 * Ku;
  parameter Real Ti = 0.5 * Tu;
  parameter Real Td = 0.33 * Tu;
  Modelica.Blocks.Continuous.LimPID pid( Td = Td, Ti = Ti, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.InitPID.NoInit,k = Kp, limitsAtInit = true, yMax = 50e-5) annotation(
    Placement(visible = true, transformation(origin = {-10, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  /*1e-12*/
  /*0.1*/
  Modelica.Blocks.Sources.Constant const(k = 10e6) annotation(
    Placement(visible = true, transformation(origin = {-70, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp ramp1(duration = 2, height = -1e6, offset = 11000000.000000002) annotation(
    Placement(visible = true, transformation(origin = {-70, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(greaterThreshold.y, switch1.u2) annotation(
    Line(points = {{-17, 50}, {28, 50}}, color = {255, 0, 255}));
  connect(ramp.y, switch1.u1) annotation(
    Line(points = {{-59, 90}, {0, 90}, {0, 58}, {28, 58}}, color = {0, 0, 127}));
  connect(step.y, greaterThreshold.u) annotation(
    Line(points = {{-59, 50}, {-40, 50}}, color = {0, 0, 127}));
  connect(pointKinetics.power, switch1.u3) annotation(
    Line(points = {{10, 20}, {20, 20}, {20, 42}, {28, 42}}, color = {0, 0, 127}));
  connect(pointKinetics.externalReactivity, externalReactivityController.externalReactivity) annotation(
    Line(points = {{-10, 20}, {-40, 20}, {-40, -20}, {-10, -20}}, color = {0, 0, 127}));
  connect(const.y, pid.u_s) annotation(
    Line(points = {{-59, -50}, {-22, -50}}, color = {0, 0, 127}));
  connect(switch1.y, externalReactivityController.power) annotation(
    Line(points = {{52, 50}, {60, 50}, {60, -20}, {12, -20}}, color = {0, 0, 127}));
  connect(ramp1.y, pid.u_m) annotation(
    Line(points = {{-58, -90}, {-10, -90}, {-10, -62}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end ExternalReactivityControllerSimulator;