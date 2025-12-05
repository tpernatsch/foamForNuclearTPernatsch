model DoublingTimeCalculator
  Modelica.Blocks.Interfaces.RealOutput doublingTime annotation(
    Placement(transformation(origin = {150, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput power annotation(
    Placement(transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Logical.Switch switch1 annotation(
    Placement(transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.GreaterEqual greaterEqual annotation(
    Placement(transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y = 1e-6) annotation(
    Placement(transformation(origin = {-90, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Division division annotation(
    Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain gain(k = log(2)) annotation(
    Placement(transformation(origin = {70, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Derivative derivative annotation(
    Placement(transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.Or or1 annotation(
    Placement(transformation(origin = {20, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.LessEqual lessEqual annotation(
    Placement(transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain gain1(k = -1)  annotation(
    Placement(transformation(origin = {-50, -38}, extent = {{-4, -4}, {4, 4}})));
equation
  connect(realExpression.y, greaterEqual.u2) annotation(
    Line(points = {{-79, -50}, {-30, -50}, {-30, -8}, {-22, -8}}, color = {0, 0, 127}, thickness = 0.5));
  connect(realExpression.y, switch1.u3) annotation(
    Line(points = {{-79, -50}, {41, -50}, {41, -8}, {58, -8}}, color = {0, 0, 127}, thickness = 0.5));
  connect(power, gain.u) annotation(
    Line(points = {{-120, 0}, {-90, 0}, {-90, 50}, {58, 50}}, color = {0, 0, 127}, thickness = 0.5));
  connect(switch1.y, division.u2) annotation(
    Line(points = {{81, 0}, {89, 0}, {89, -6}, {97, -6}}, color = {0, 0, 127}, thickness = 0.5));
  connect(gain.y, division.u1) annotation(
    Line(points = {{81, 50}, {89, 50}, {89, 6}, {97, 6}}, color = {0, 0, 127}, thickness = 0.5));
  connect(division.y, doublingTime) annotation(
    Line(points = {{121, 0}, {149, 0}}, color = {0, 0, 127}, thickness = 0.5));
  connect(power, derivative.u) annotation(
    Line(points = {{-120, 0}, {-82, 0}}, color = {0, 0, 127}));
  connect(derivative.y, greaterEqual.u1) annotation(
    Line(points = {{-59, 0}, {-22, 0}}, color = {0, 0, 127}));
  connect(derivative.y, switch1.u1) annotation(
    Line(points = {{-59, 0}, {-30, 0}, {-30, 20}, {40, 20}, {40, 8}, {58, 8}}, color = {0, 0, 127}));
  connect(or1.y, switch1.u2) annotation(
    Line(points = {{32, 0}, {58, 0}}, color = {255, 0, 255}));
  connect(greaterEqual.y, or1.u1) annotation(
    Line(points = {{2, 0}, {8, 0}}, color = {255, 0, 255}));
  connect(lessEqual.y, or1.u2) annotation(
    Line(points = {{2, -30}, {8, -30}, {8, -8}}, color = {255, 0, 255}));
  connect(derivative.y, lessEqual.u1) annotation(
    Line(points = {{-58, 0}, {-40, 0}, {-40, -30}, {-22, -30}}, color = {0, 0, 127}));
  connect(gain1.y, lessEqual.u2) annotation(
    Line(points = {{-46, -38}, {-22, -38}}, color = {0, 0, 127}));
  connect(realExpression.y, gain1.u) annotation(
    Line(points = {{-78, -50}, {-68, -50}, {-68, -38}, {-55, -38}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")));
end DoublingTimeCalculator;