model SignalProcess
  parameter Real powerGain = 1;
  parameter Real minValue = 0;
  parameter Real maxValue = 1;
  Modelica.Blocks.Math.Gain gain1(k = powerGain) annotation(
    Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput in_power annotation(
    Placement(visible = true, transformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Abs abs1 annotation(
    Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput out_power annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(limitsAtInit = true, uMax = maxValue, uMin = minValue)  annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  /*
  Modelica.Blocks.Nonlinear.SlewRateLimiter slewRateLimiter(Rising = 1000)  annotation(
    Placement(visible = true, transformation(origin = {-36, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  */
equation
  connect(gain1.y, out_power) annotation(
    Line(points = {{82, 0}, {110, 0}}, color = {0, 0, 127}));
  connect(in_power, abs1.u) annotation(
    Line(points = {{-110, 0}, {-82, 0}}, color = {0, 0, 127}));
  connect(limiter.y, gain1.u) annotation(
    Line(points = {{12, 0}, {58, 0}}, color = {0, 0, 127}));
  connect(abs1.y, limiter.u) annotation(
    Line(points = {{-58, 0}, {-12, 0}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end SignalProcess;