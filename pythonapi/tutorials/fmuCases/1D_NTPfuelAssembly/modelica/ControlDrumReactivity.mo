model ControlDrumReactivity
  Modelica.Blocks.Math.Tanh tanh annotation(
    Placement(transformation(origin = {-8, -6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain gain(k = 0.020) annotation(
    Placement(transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add annotation(
    Placement(transformation(origin = {-38, -6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y = -1.630) annotation(
    Placement(transformation(origin = {-78, -24}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add1 annotation(
    Placement(transformation(origin = {62, -12}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(y = -1.083) annotation(
    Placement(transformation(origin = {22, -40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain gain1(k = 3.228) annotation(
    Placement(transformation(origin = {22, -6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain convertDollarToReact(k = 2*653.51e-5) annotation(
    Placement(transformation(origin = {102, -12}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput angle annotation(
    Placement(transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput reactivity annotation(
    Placement(transformation(origin = {140, -12}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(gain.y, add.u1) annotation(
    Line(points = {{-67, 0}, {-51, 0}}, color = {0, 0, 127}));
  connect(realExpression.y, add.u2) annotation(
    Line(points = {{-67, -24}, {-61, -24}, {-61, -12}, {-51, -12}}, color = {0, 0, 127}));
  connect(add.y, tanh.u) annotation(
    Line(points = {{-27, -6}, {-21, -6}}, color = {0, 0, 127}));
  connect(gain1.y, add1.u1) annotation(
    Line(points = {{33, -6}, {50, -6}}, color = {0, 0, 127}));
  connect(gain1.u, tanh.y) annotation(
    Line(points = {{10, -6}, {3, -6}}, color = {0, 0, 127}));
  connect(realExpression1.y, add1.u2) annotation(
    Line(points = {{33, -40}, {42, -40}, {42, -18}, {50, -18}}, color = {0, 0, 127}));
  connect(add1.y, convertDollarToReact.u) annotation(
    Line(points = {{73, -12}, {89, -12}}, color = {0, 0, 127}));
  connect(angle, gain.u) annotation(
    Line(points = {{-120, 0}, {-90, 0}}, color = {0, 0, 127}));
  connect(reactivity, convertDollarToReact.y) annotation(
    Line(points = {{140, -12}, {114, -12}}, color = {0, 0, 127}));

annotation(
    uses(Modelica(version = "3.2.3")),
  Documentation(info = "<html><head></head><body>Reactivity function of the control drum with \"angle\" in deg and \"reactivity\" in unit.</body></html>"),
  Diagram(coordinateSystem(extent = {{-140, 20}, {160, -60}})),
  version = "");
end ControlDrumReactivity;