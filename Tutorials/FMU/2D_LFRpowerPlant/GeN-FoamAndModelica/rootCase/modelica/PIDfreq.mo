model PIDfreq
  parameter Real Ku = 2;
  parameter Real Tu = 0.5;
  parameter Real Kp = 0.2*Ku;
  parameter Real Ti = 0.5*Tu;
  parameter Real Td = 0.33*Tu*1000;
  Modelica.Blocks.Continuous.LimPID pidFrequency(
    k = Kp, Td = Td, Ti = Ti, controllerType = Modelica.Blocks.Types.SimpleController.PI
  ) annotation(
    Placement(visible = true, transformation(origin = {0, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput measFrequency annotation(
    Placement(visible = true, transformation(origin = {-120, 30}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput pidFrequencyResponse annotation(
    Placement(visible = true, transformation(origin = {110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = 50)  annotation(
    Placement(visible = true, transformation(origin = {-70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //
  //
  // Ziegler–Nichols tuning method (no overshoot)
  parameter Real KuReact = 1e-12;
  parameter Real TuReact = 5/2;
  parameter Real KpReact = 0.2*KuReact;
  parameter Real TiReact = 0.5*TuReact;
  parameter Real TdReact = 0.33*TuReact;
  Modelica.Blocks.Continuous.LimPID pidExternalSource(Td = TdExtSrc, Ti = TiExtSrc, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.InitPID.InitialOutput, k = KpExtSrc, limitsAtInit = true, yMax = 1.2, yMin = 0, y_start = 1) annotation(
    Placement(visible = true, transformation(origin = {60, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add addPowers annotation(
    Placement(visible = true, transformation(origin = {-140, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput externalSourceModulation(start = 1) annotation(
    Placement(visible = true, transformation(origin = {110, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain machineEfficiency(k = 300 / 120) annotation(
    Placement(visible = true, transformation(origin = {-110, -130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput outerCorePower(start = 13755) annotation(
    Placement(visible = true, transformation(origin = {-180, -100}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Math.Gain corePowerMeasure(k = 72 * 198.28145) annotation(
    Placement(visible = true, transformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput innerCorePower(start = 6487) annotation(
    Placement(visible = true, transformation(origin = {-180, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput measElectricalPower(start = 6487) annotation(
    Placement(visible = true, transformation(origin = {-180, -130}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput externalReactivity(start = 0) annotation(
    Placement(visible = true, transformation(origin = {110, -130}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //
  // Ziegler–Nichols tuning method (no overshoot)
  parameter Real KuExtSrc = 1e-9;
  parameter Real TuExtSrc = 5/2;
  parameter Real KpExtSrc = 0.2*KuExtSrc;
  parameter Real TiExtSrc = 0.5*TuExtSrc;
  parameter Real TdExtSrc = 0.33*TuExtSrc;
  Modelica.Blocks.Continuous.LimPID pidReactivity(Td = TdReact, Ti = TiReact, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.InitPID.InitialOutput, k = KpReact, limitsAtInit = true, yMax = 250e-5, y_start = 0) annotation(
    Placement(visible = true, transformation(origin = {60, -130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(measFrequency, pidFrequency.u_m) annotation(
    Line(points = {{-120, 30}, {0, 30}, {0, 58}}, color = {0, 0, 127}));
  connect(const.y, pidFrequency.u_s) annotation(
    Line(points = {{-59, 70}, {-12, 70}}, color = {0, 0, 127}));
  connect(pidFrequency.y, pidFrequencyResponse) annotation(
    Line(points = {{11, 70}, {110, 70}}, color = {0, 0, 127}));
  connect(addPowers.u1, innerCorePower) annotation(
    Line(points = {{-152, -74}, {-161, -74}, {-161, -60}, {-180, -60}}, color = {0, 0, 127}));
  connect(pidExternalSource.y, externalSourceModulation) annotation(
    Line(points = {{71, -10}, {110, -10}}, color = {0, 0, 127}));
  connect(machineEfficiency.y, pidExternalSource.u_s) annotation(
    Line(points = {{-99, -130}, {-20, -130}, {-20, -10}, {48, -10}}, color = {204, 0, 0}));
  connect(addPowers.y, corePowerMeasure.u) annotation(
    Line(points = {{-129, -80}, {-122, -80}}, color = {0, 0, 127}));
  connect(addPowers.u2, outerCorePower) annotation(
    Line(points = {{-152, -86}, {-161, -86}, {-161, -100}, {-180, -100}}, color = {0, 0, 127}));
  connect(measElectricalPower, machineEfficiency.u) annotation(
    Line(points = {{-180, -130}, {-122, -130}}, color = {0, 0, 127}));
  connect(pidReactivity.y, externalReactivity) annotation(
    Line(points = {{71, -130}, {110, -130}}, color = {0, 0, 127}));
  connect(machineEfficiency.y, pidReactivity.u_s) annotation(
    Line(points = {{-99, -130}, {48, -130}}, color = {255, 0, 0}));
  connect(corePowerMeasure.y, pidExternalSource.u_m) annotation(
    Line(points = {{-98, -80}, {60, -80}, {60, -22}}, color = {0, 0, 127}));
  connect(corePowerMeasure.y, pidReactivity.u_m) annotation(
    Line(points = {{-98, -80}, {20, -80}, {20, -160}, {60, -160}, {60, -142}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "3.2.3")),
  Diagram(coordinateSystem(extent = {{-200, 80}, {120, -180}})),
  version = "");
end PIDfreq;