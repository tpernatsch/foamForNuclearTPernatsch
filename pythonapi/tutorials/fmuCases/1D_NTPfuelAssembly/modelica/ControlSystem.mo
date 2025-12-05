model ControlSystem
  parameter Real nElement = 1500;
  parameter Real wnom;
  parameter Real burnTimeFullThrottle = 620;
  Modelica.Blocks.Sources.Ramp powerCmd(duration = 10, height = -737e6, offset = 937e6, startTime = 20) annotation(
    Placement(transformation(origin = {-80, -120}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput extReactivity annotation(
    Placement(transformation(origin = {280, -90}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {210, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput power annotation(
    Placement(transformation(origin = {-120, -160}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-220, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Math.Gain corePower(k = nElement) annotation(
    Placement(transformation(origin = {-40, -160}, extent = {{-10, -10}, {10, 10}})));
  ControlDrumReactivity controlDrumReactivity annotation(
    Placement(transformation(origin = {240, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add annotation(
    Placement(transformation(origin = {150, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression offsetCD(y = 99.9) annotation(
    Placement(transformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Nonlinear.Limiter limiterCDangle(uMax = 180, uMin = 0) annotation(
    Placement(transformation(origin = {180, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput massFlow annotation(
    Placement(transformation(origin = {-120, 30}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-220, 100}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Sources.Ramp massFlowCmd(height = 0, duration = 100, offset = wnom, startTime = 1000000) annotation(
    Placement(transformation(origin = {-130, 62}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput turbineValveOpening annotation(
    Placement(transformation(origin = {10, 70}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {210, 100}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput Tnozzle annotation(
    Placement(transformation(origin = {-220, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-220, -100}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Logical.Switch switchPowerCmd annotation(
    Placement(transformation(origin = {-30, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.BooleanExpression selectPowerFromThrust(y = true) annotation(
    Placement(transformation(origin = {-90, -90}, extent = {{-10, -10}, {10, 10}})));
  ThrustCalculator thrustCalculator annotation(
    Placement(transformation(origin = {-10, 6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp thrustStartRamp(height = 220e3, duration = 20, offset = 10e3, startTime = 10) annotation(
    Placement(transformation(origin = {-250, -100}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp thrustEndRamp(height = -220e3, duration = 20, offset = 0, startTime = 10 + 20 + burnTimeFullThrottle) annotation(
    Placement(transformation(origin = {-250, -140}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add thrustCmd annotation(
    Placement(transformation(origin = {-210, -120}, extent = {{-10, -10}, {10, 10}})));
  DoublingTimeCalculator doublingTimeCalculator annotation(
    Placement(transformation(origin = {0, -160}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID PI_Thrust(controllerType = Modelica.Blocks.Types.SimpleController.PI, k = 30, Ti = 0.1, yMax = 500e6, yMin = 100) annotation(
    Placement(transformation(origin = {-90, -60}, extent = {{-10, 10}, {10, -10}})));
  DeltaVCalculator deltaVCalculator(propellantMass = 34602.7)  annotation(
    Placement(transformation(origin = {70, 6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID PI_massFlow(controllerType = Modelica.Blocks.Types.SimpleController.PI, Ti = 0.1, k = -0.01, initType = Modelica.Blocks.Types.InitPID.InitialOutput, y_start = 0.5, yMax = 1, yMin = 0.01)  annotation(
    Placement(transformation(origin = {-40, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression minDoublingTime(y = 20)  annotation(
    Placement(transformation(origin = {0, -140}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID PI_DoublingTime(controllerType = Modelica.Blocks.Types.SimpleController.PI, Ti = 0.1, yMax = 80, yMin = -100, k = -0.01)  annotation(
    Placement(transformation(origin = {50, -140}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID PI_Power(controllerType = Modelica.Blocks.Types.SimpleController.PI, k = 5e-9, Ti = 0.05, yMax = 80, yMin = -100, Td = 1000)  annotation(
    Placement(transformation(origin = {10, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Min min annotation(
    Placement(transformation(origin = {110, -96}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrderTd(T = 1)  annotation(
    Placement(transformation(origin = {40, -160}, extent = {{-4, -4}, {4, 4}})));
  Modelica.Blocks.Logical.Switch switchThrustCmd annotation(
    Placement(transformation(origin = {-150, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.BooleanExpression selectThrustFromDeltaV(y = true) annotation(
    Placement(transformation(origin = {-190, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID PI_DeltaV(Ti = 0.1, controllerType = Modelica.Blocks.Types.SimpleController.PI, k = 100, yMax = 239e3, yMin = 0) annotation(
    Placement(transformation(origin = {-210, -60}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Blocks.Sources.Step deltaVcmd(height = 3119, startTime = 0)  annotation(
    Placement(transformation(origin = {-270, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Nonlinear.Limiter limiterTd(uMax = 1e9, uMin = 19)  annotation(
    Placement(transformation(origin = {20, -160}, extent = {{-4, -4}, {4, 4}})));
  Modelica.Blocks.Nonlinear.SlewRateLimiter slewRateLimiterCDangle(Rising = 0.5, Falling = -5)  annotation(
    Placement(transformation(origin = {210, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.BooleanExpression selectMassFlowFromTnozzle(y = true) annotation(
    Placement(transformation(origin = {-100, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.Switch switch annotation(
    Placement(transformation(origin = {-70, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.LimPID PI_Tnozzle(Ti = 1, controllerType = Modelica.Blocks.Types.SimpleController.PI, initType = Modelica.Blocks.Types.InitPID.InitialOutput, k = -0.001, yMax = wnom, yMin = 0.5*wnom, y_start = wnom) annotation(
    Placement(transformation(origin = {-190, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression TnozzleCmd(y = 500)  annotation(
    Placement(transformation(origin = {-230, 90}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(offsetCD.y, add.u1) annotation(
    Line(points = {{121, -60}, {129, -60}, {129, -84}, {137, -84}}, color = {0, 0, 127}, thickness = 0.5));
  connect(add.y, limiterCDangle.u) annotation(
    Line(points = {{161, -90}, {168, -90}}, color = {0, 0, 127}, thickness = 0.5));
  connect(powerCmd.y, switchPowerCmd.u3) annotation(
    Line(points = {{-69, -120}, {-60, -120}, {-60, -98}, {-42, -98}}, color = {0, 0, 127}, thickness = 0.5));
  connect(selectPowerFromThrust.y, switchPowerCmd.u2) annotation(
    Line(points = {{-79, -90}, {-42, -90}}, color = {255, 0, 255}, thickness = 0.5));
  connect(Tnozzle, thrustCalculator.Tnozzle) annotation(
    Line(points = {{-220, 0}, {-22, 0}}, color = {0, 0, 127}, thickness = 0.5));
  connect(massFlow, thrustCalculator.massFlow) annotation(
    Line(points = {{-120, 30}, {-40, 30}, {-40, 12}, {-22, 12}}, color = {0, 0, 127}));
  connect(thrustCalculator.thrust, PI_Thrust.u_m) annotation(
    Line(points = {{1, 0}, {19, 0}, {19, -22}, {-90, -22}, {-90, -48}}, color = {0, 0, 127}, thickness = 0.5));
  connect(thrustStartRamp.y, thrustCmd.u1) annotation(
    Line(points = {{-239, -100}, {-232, -100}, {-232, -114}, {-222, -114}}, color = {0, 0, 127}, thickness = 0.5));
  connect(thrustEndRamp.y, thrustCmd.u2) annotation(
    Line(points = {{-239, -140}, {-232, -140}, {-232, -126}, {-222, -126}}, color = {0, 0, 127}, thickness = 0.5));
  connect(PI_Thrust.y, switchPowerCmd.u1) annotation(
    Line(points = {{-79, -60}, {-60, -60}, {-60, -82}, {-42, -82}}, color = {0, 0, 127}, thickness = 0.5));
  connect(thrustCalculator.thrust, deltaVCalculator.thrust) annotation(
    Line(points = {{2, 0}, {58, 0}}, color = {0, 0, 127}));
  connect(massFlow, deltaVCalculator.massFlow) annotation(
    Line(points = {{-120, 30}, {42, 30}, {42, 10}, {58, 10}}, color = {0, 0, 127}));
  connect(PI_massFlow.u_m, massFlow) annotation(
    Line(points = {{-40, 58}, {-40, 30}, {-120, 30}}, color = {0, 0, 127}));
  connect(PI_massFlow.y, turbineValveOpening) annotation(
    Line(points = {{-29, 70}, {10, 70}}, color = {0, 0, 127}));
  connect(corePower.y, doublingTimeCalculator.power) annotation(
    Line(points = {{-29, -160}, {-12, -160}}, color = {0, 0, 127}, thickness = 0.5));
  connect(switchPowerCmd.y, PI_Power.u_s) annotation(
    Line(points = {{-18, -90}, {-2, -90}}, color = {0, 0, 127}));
  connect(corePower.y, PI_Power.u_m) annotation(
    Line(points = {{-28, -160}, {-20, -160}, {-20, -110}, {10, -110}, {10, -102}}, color = {0, 0, 127}, thickness = 0.5));
  connect(minDoublingTime.y, PI_DoublingTime.u_s) annotation(
    Line(points = {{12, -140}, {38, -140}}, color = {0, 0, 127}, thickness = 0.5));
  connect(PI_Power.y, min.u1) annotation(
    Line(points = {{22, -90}, {98, -90}}, color = {0, 0, 127}, thickness = 0.5));
  connect(min.y, add.u2) annotation(
    Line(points = {{122, -96}, {138, -96}}, color = {0, 0, 127}, thickness = 0.5));
  connect(PI_DoublingTime.y, min.u2) annotation(
    Line(points = {{62, -140}, {80, -140}, {80, -102}, {98, -102}}, color = {0, 0, 127}, thickness = 0.5));
  connect(selectThrustFromDeltaV.y, switchThrustCmd.u2) annotation(
    Line(points = {{-179, -90}, {-163, -90}}, color = {255, 0, 255}));
  connect(thrustCmd.y, switchThrustCmd.u3) annotation(
    Line(points = {{-199, -120}, {-180, -120}, {-180, -98}, {-162, -98}}, color = {0, 0, 127}, thickness = 0.5));
  connect(PI_DeltaV.y, switchThrustCmd.u1) annotation(
    Line(points = {{-198, -60}, {-170, -60}, {-170, -82}, {-162, -82}}, color = {0, 0, 127}, thickness = 0.5));
  connect(switchThrustCmd.y, PI_Thrust.u_s) annotation(
    Line(points = {{-138, -90}, {-120, -90}, {-120, -60}, {-102, -60}}, color = {0, 0, 127}, thickness = 0.5));
  connect(deltaVCalculator.deltaV, PI_DeltaV.u_m) annotation(
    Line(points = {{82, 6}, {100, 6}, {100, -12}, {-210, -12}, {-210, -48}}, color = {0, 0, 127}, thickness = 0.5));
  connect(deltaVcmd.y, PI_DeltaV.u_s) annotation(
    Line(points = {{-258, -60}, {-222, -60}}, color = {0, 0, 127}, thickness = 0.5));
  connect(controlDrumReactivity.reactivity, extReactivity) annotation(
    Line(points = {{251, -90}, {279, -90}}, color = {0, 0, 127}));
  connect(power, corePower.u) annotation(
    Line(points = {{-120, -160}, {-52, -160}}, color = {0, 0, 127}));
  connect(doublingTimeCalculator.doublingTime, limiterTd.u) annotation(
    Line(points = {{12, -160}, {16, -160}}, color = {0, 0, 127}));
  connect(firstOrderTd.y, PI_DoublingTime.u_m) annotation(
    Line(points = {{44, -160}, {50, -160}, {50, -152}}, color = {0, 0, 127}));
  connect(limiterTd.y, firstOrderTd.u) annotation(
    Line(points = {{24, -160}, {35, -160}}, color = {0, 0, 127}));
  connect(limiterCDangle.y, slewRateLimiterCDangle.u) annotation(
    Line(points = {{192, -90}, {198, -90}}, color = {0, 0, 127}));
  connect(slewRateLimiterCDangle.y, controlDrumReactivity.angle) annotation(
    Line(points = {{222, -90}, {228, -90}}, color = {0, 0, 127}));
  connect(selectMassFlowFromTnozzle.y, switch.u2) annotation(
    Line(points = {{-89, 70}, {-82, 70}}, color = {255, 0, 255}));
  connect(massFlowCmd.y, switch.u3) annotation(
    Line(points = {{-119, 62}, {-82, 62}}, color = {0, 0, 127}));
  connect(switch.y, PI_massFlow.u_s) annotation(
    Line(points = {{-58, 70}, {-52, 70}}, color = {0, 0, 127}));
  connect(PI_Tnozzle.y, switch.u1) annotation(
    Line(points = {{-178, 90}, {-100, 90}, {-100, 78}, {-82, 78}}, color = {0, 0, 127}));
  connect(Tnozzle, PI_Tnozzle.u_m) annotation(
    Line(points = {{-220, 0}, {-190, 0}, {-190, 78}}, color = {0, 0, 127}));
  connect(TnozzleCmd.y, PI_Tnozzle.u_s) annotation(
    Line(points = {{-218, 90}, {-202, 90}}, color = {0, 0, 127}));
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -200}, {200, 200}}, initialScale = 0.1), graphics = {Rectangle(lineColor = {170, 170, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-200, 200}, {200, -200}}), Text(textColor = {170, 170, 255}, extent = {{-140, 140}, {140, -140}}, textString = "C")}),
    uses(Modelica(version = "3.2.3")),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection --allowNonStandardModelica=reinitInAlgorithms -d=initialization,NLSanalyticJacobian",
    Diagram(coordinateSystem(extent = {{-280, 100}, {300, -180}})),
    version = "",
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
end ControlSystem;