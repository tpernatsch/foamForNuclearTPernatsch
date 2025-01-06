model SecondaryCircuitWithSimulink
  // package Medium = Modelica.Media.Water.WaterIF97_ph;
  replaceable package Water = ThermoPower.Water.StandardWater constrainedby Modelica.Media.Interfaces.PartialPureSubstance "Fluid model";
  constant Real pi = Modelica.Constants.pi;
  parameter Real wnom = 192.7;
  constant Real bypassFrac = 0;
  parameter Modelica.SIunits.Pressure maxP = 18.8e6 "Maximum pressure";
  parameter Modelica.SIunits.Pressure interP = 6e6 "Intermediate pressure";
  parameter Modelica.SIunits.Pressure minP = 0.008e6 "Minimum pressure";
  parameter Integer nSGparallel = 8;
  parameter Integer nBayonnets_ = 510*nSGparallel;
  parameter Real pinPitch = 31.73e-3*1.42;
  parameter Real outerTubeDiameter = 25.40e-3;
  parameter Real innerTubeDiameter = 19.05e-3;
  parameter Real volumetricArea = 56.69862 "m2/m3";
  // 48.94457;
  parameter Real Pm_ = pi*(outerTubeDiameter + innerTubeDiameter) "Wet perimeter";
  parameter Real Dh_ = outerTubeDiameter - innerTubeDiameter "Hydraulic diameter for a ring";
  parameter Real A_ = sqrt(3.0)*pinPitch/4.0 "Cross-sectional area for a single pin (pattern of the assembly)";
  parameter Real Pc_ = pi*outerTubeDiameter "Heated perimeter";
  parameter Real powerGain = 360.0/5.0*volumetricArea;
  inner ThermoPower.System system annotation(
    Placement(visible = true, transformation(origin = {50, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  /*pi*(innerTubeDiameter^2)/4.0*/
  ThermoPower.Water.FlowSplit flowSplitTurbLP annotation(
    Placement(visible = true, transformation(origin = {30, 100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  ThermoPower.Water.ValveVap valveVapAdmissionTurb(Av = 0.2, CheckValve = true, dpnom = 999.9999999999999, pnom = maxP, rhonom = 0.08e3, theta_fix = 1, useThetaInput = true, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-90, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Water.SteamTurbineStodola HP(Kt = 0.0055, PRstart = maxP/interP, pin(fixed = false), pnom = maxP, pout(fixed = false), wnom = wnom, wstart = wnom, eta_mech = 0.7, eta_iso_nom = 0.7) annotation(
    Placement(visible = true, transformation(origin = {-50, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  ThermoPower.Water.SteamTurbineStodola LP1(Kt = 0.00929, PRstart = interP/minP, pin(fixed = false), pnom = interP, pout(fixed = false), wnom = wnom/2, eta_mech = 0.7, eta_iso_nom = 0.7) annotation(
    Placement(visible = true, transformation(origin = {0, 60}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
  ThermoPower.Water.SteamTurbineStodola LP2(Kt = 0.00929, PRstart = interP/minP, pin(fixed = false), pnom = interP, pout(fixed = false), wnom = wnom/2, eta_mech = 0.7, eta_iso_nom = 0.7) annotation(
    Placement(visible = true, transformation(origin = {60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  ThermoPower.Water.FlowJoin flowJoinTurbLP annotation(
    Placement(visible = true, transformation(origin = {30, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  ThermoPower.Water.SinkPressure sinkPressure(T = 41 + 273, p0 = minP, use_T = true) annotation(
    Placement(visible = true, transformation(origin = {90, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Electrical.Generator generator(Pnom(displayUnit = "MW") = 125e6, referenceGenerator = true) annotation(
    Placement(visible = true, transformation(origin = {160, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Electrical.Load load(Pnom(displayUnit = "MW") = 125e6, usePowerInput = true) annotation(
    Placement(visible = true, transformation(origin = {190, -130}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  ThermoPower.Electrical.PowerSensor powerSensorElectric annotation(
    Placement(visible = true, transformation(origin = {190, -50}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
  Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensorMech annotation(
    Placement(visible = true, transformation(origin = {130, 60}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  ThermoPower.Water.PressDropLin pressDropLin(R = 1/wnom) annotation(
    Placement(visible = true, transformation(origin = {-10, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Steam generator
  parameter Real wnm = 1;
  ThermoPower.Water.Flow1DFV steamGenBayonnet(A = A_, Dhyd = innerTubeDiameter, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.TwoPhases, H = -6, L = 6, Nt = nBayonnets_, hstartin = 1e6, hstartout = 1.540e6, noInitialPressure = true, omega = 0, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-100, -90}, extent = {{10, 10}, {-10, -10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG1(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid, H = 0.82, L = 0.82, Nt = nBayonnets_, hstartin = 1.541e6, hstartout = 1.732e6, noInitialPressure = true, omega = Pc_, p(fixed = false), pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, -170}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG2(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Liquid, H = 0.59, L = 0.59, Nt = nBayonnets_, hstartin = 1.732e6, hstartout = 1.889e6, noInitialPressure = true, omega = Pc_, p(fixed = false), pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG3(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.TwoPhases, H = 0.3, L = 0.3, Nt = nBayonnets_, hstartin = 1.889e6, hstartout = 1.978e6, noInitialPressure = true, omega = Pc_, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG4(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.TwoPhases, H = 0.4, L = 0.4, Nt = nBayonnets_, hstartin = 1.978e6, hstartout = 2.106e6, noInitialPressure = true, omega = Pc_, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG5(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.TwoPhases, H = 0.69, L = 0.69, Nt = nBayonnets_, hstartin = 2.106e6, hstartout = 2.361e6, noInitialPressure = true, omega = Pc_, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG6(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.TwoPhases, H = 0.90, L = 0.90, Nt = nBayonnets_, hstartin = 2.361e6, hstartout = 2.689e6, noInitialPressure = true, omega = Pc_, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG7(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Steam, H = 1.5, L = 1.5, Nt = nBayonnets_, hstartin = 2.689e6, hstartout = 2.985e6, noInitialPressure = true, omega = Pc_, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Water.Flow1DFV2ph SG8(A = A_, Dhyd = Dh_, FluidPhaseStart = ThermoPower.Choices.FluidPhase.FluidPhases.Steam, H = 0.80, L = 0.80, Nt = nBayonnets_, hstartin = 2.985e6, hstartout = 3.083e6, noInitialPressure = true, omega = Pc_, pstart = maxP, wnm = wnm, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-180, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  // Thermal boundaries
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG1 annotation(
    Placement(visible = true, transformation(origin = {-200, -170}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG2 annotation(
    Placement(visible = true, transformation(origin = {-200, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG3 annotation(
    Placement(visible = true, transformation(origin = {-200, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG4 annotation(
    Placement(visible = true, transformation(origin = {-200, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG5 annotation(
    Placement(visible = true, transformation(origin = {-200, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG6 annotation(
    Placement(visible = true, transformation(origin = {-200, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG7 annotation(
    Placement(visible = true, transformation(origin = {-200, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  ThermoPower.Thermal.HeatSource1DFV heatSourceSG8 annotation(
    Placement(visible = true, transformation(origin = {-200, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  // Power inputs
  Modelica.Blocks.Interfaces.RealInput steamGenPower1(start = 9080) annotation(
    Placement(visible = true, transformation(origin = {-260, -170}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, -190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower2(start = 7510) annotation(
    Placement(visible = true, transformation(origin = {-260, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, -150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower3(start = 4230) annotation(
    Placement(visible = true, transformation(origin = {-260, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower4(start = 6110) annotation(
    Placement(visible = true, transformation(origin = {-260, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower5(start = 12130) annotation(
    Placement(visible = true, transformation(origin = {-260, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower6(start = 15640) annotation(
    Placement(visible = true, transformation(origin = {-260, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower7(start = 14100) annotation(
    Placement(visible = true, transformation(origin = {-260, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput steamGenPower8(start = 4660) annotation(
    Placement(visible = true, transformation(origin = {-260, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-210, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Temperature outputs
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature1(start = 631) annotation(
    Placement(visible = true, transformation(origin = {-150, -170}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature2(start = 634) annotation(
    Placement(visible = true, transformation(origin = {-150, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature3(start = 634) annotation(
    Placement(visible = true, transformation(origin = {-150, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature4(start = 634) annotation(
    Placement(visible = true, transformation(origin = {-150, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature5(start = 634) annotation(
    Placement(visible = true, transformation(origin = {-150, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature6(start = 649) annotation(
    Placement(visible = true, transformation(origin = {-150, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature7(start = 699) annotation(
    Placement(visible = true, transformation(origin = {-150, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput steamGenTemperature8(start = 722) annotation(
    Placement(visible = true, transformation(origin = {-150, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Signal proccessing
  SignalProcess signalProcess1(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, -170}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess2(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, -140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess3(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess4(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess5(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess6(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess7(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SignalProcess signalProcess8(maxValue = 30e6, powerGain = powerGain) annotation(
    Placement(visible = true, transformation(origin = {-226, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Eletric power output
  Modelica.SIunits.Power injectedPower;
  Modelica.Blocks.Sources.Ramp rampLoad(duration = 10, height = -20e6, offset = 0, startTime = 1e15) annotation(
    Placement(visible = true, transformation(origin = {120, -130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Electrical.FrequencySensor frequencySensor annotation(
    Placement(visible = true, transformation(origin = {150, 150}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  ThermoPower.Water.SourcePressure sourcePressure(T = 335 + 273.15, p0(displayUnit = "Pa") = maxP, use_T = true, use_in_p0 = false) annotation(
    Placement(visible = true, transformation(origin = {-30, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  ThermoPower.Water.SensW sensW annotation(
    Placement(visible = true, transformation(origin = {-164, 68}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  ThermoPower.Water.FlowJoin flowJoin annotation(
    Placement(visible = true, transformation(origin = {50, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Water.FlowSplit flowSplitBypass annotation(
    Placement(visible = true, transformation(origin = {-120, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Water.ValveVap valveVapBypass(Av = 0.5, CheckValve = true, dpnom = 10000, pnom = maxP, rhonom = 0.08e3, theta_fix = 0, useThetaInput = false, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-30, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Water.ValveLiqChoked valveLiqChoked(Av = 0.1, CheckValve = true, dpnom = 100, pnom = maxP, theta_fix = 1, useThetaInput = false, wnom = wnom) annotation(
    Placement(visible = true, transformation(origin = {-70, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 1000) annotation(
    Placement(visible = true, transformation(origin = {98, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //
  // Ziegler–Nichols tuning method (no overshoot)
  parameter Real Ku = 0.25;
  // /5;
  parameter Real Tu = 1/10*2;
  // *10;
  parameter Real Kp = 0.2*Ku;
  // 0.6*Ku; //
  parameter Real Ti = 0.5*Tu;
  // 0.5*Tu; //
  parameter Real Td = 0.33*Tu*1000;
  // 0.125*Tu; //
  /*10000*/
  /*0.1*/
  /*0.1*/
  //
  // Ziegler–Nichols tuning method (no overshoot)
  parameter Real KuReact = 1e-12;
  parameter Real TuReact = 5/2;
  parameter Real KpReact = 0.2*KuReact;
  parameter Real TiReact = 0.5*TuReact;
  parameter Real TdReact = 0.33*TuReact;
  Modelica.Blocks.Continuous.LimPID pidReactivity(Td = TdReact, Ti = TiReact, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.InitPID.InitialOutput, k = KpReact, limitsAtInit = true, yMax = 250e-5, y_start = 0) annotation(
    Placement(visible = true, transformation(origin = {390, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step stepActivateExtReactivity(height = -1, offset = 1, startTime = 1e15) annotation(
    Placement(visible = true, transformation(origin = {250, -150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch1 annotation(
    Placement(visible = true, transformation(origin = {350, -150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = 0.5) annotation(
    Placement(visible = true, transformation(origin = {290, -150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput externalReactivity(start = 0) annotation(
    Placement(visible = true, transformation(origin = {430, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {210, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput innerCorePower(start = 6487) annotation(
    Placement(visible = true, transformation(origin = {220, 20}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-220, 140}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput outerCorePower(start = 13755) annotation(
    Placement(visible = true, transformation(origin = {220, -20}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-220, 180}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 300/120) annotation(
    Placement(visible = true, transformation(origin = {260, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add addPowers annotation(
    Placement(visible = true, transformation(origin = {260, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain corePowerMeasure(k = 72*198.28145) annotation(
    Placement(visible = true, transformation(origin = {290, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Ramp rampLoad2(duration = 10, height = 20e6, offset = 0, startTime = 1e15) annotation(
    Placement(visible = true, transformation(origin = {120, -160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add3 add3 annotation(
    Placement(visible = true, transformation(origin = {160, -130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const(k = 120e6) annotation(
    Placement(visible = true, transformation(origin = {120, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  /*
          ThermoPower.Water.SourceMassFlow sourceMassFlow(T = 335 + 273.15,p0 = maxP, use_T = true, w0 = wnom)  annotation(
            Placement(visible = true, transformation(origin = {-18, -76}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          */
  //
  //
  Modelica.Blocks.Interfaces.RealOutput externalSourceModulation(start = 1) annotation(
    Placement(visible = true, transformation(origin = {430, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {210, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //
  // Ziegler–Nichols tuning method (no overshoot)
  parameter Real KuExtSrc = 1e-9;
  parameter Real TuExtSrc = 5/2;
  parameter Real KpExtSrc = 0.2*KuExtSrc;
  parameter Real TiExtSrc = 0.5*TuExtSrc;
  parameter Real TdExtSrc = 0.33*TuExtSrc;
  Modelica.Blocks.Continuous.LimPID pidExternalSource(Td = TdExtSrc, Ti = TiExtSrc, controllerType = Modelica.Blocks.Types.SimpleController.PID, initType = Modelica.Blocks.Types.InitPID.InitialOutput, k = KpExtSrc, limitsAtInit = true, yMax = 1.2, yMin = 0, y_start = 1) annotation(
    Placement(visible = true, transformation(origin = {390, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold1(threshold = 0.5) annotation(
    Placement(visible = true, transformation(origin = {290, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step stepActivateExtSource(height = -1, offset = 1, startTime = 1e15) annotation(
    Placement(visible = true, transformation(origin = {250, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch annotation(
    Placement(visible = true, transformation(origin = {350, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput valveControlInput annotation(
    Placement(visible = true, transformation(origin = {130, 170}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {-100, 138}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput frequencyOutput(start = 50) annotation(
    Placement(visible = true, transformation(origin = {120, 150}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {50, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = 1, uMin = 0.01)  annotation(
    Placement(visible = true, transformation(origin = {50, 180}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add add annotation(
    Placement(visible = true, transformation(origin = {90, 180}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant defaultValve(k = 0.285)  annotation(
    Placement(visible = true, transformation(origin = {130, 190}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
equation
  connect(flowSplitTurbLP.out1, LP2.inlet) annotation(
    Line(points = {{34, 94}, {34, 76}, {44, 76}}, color = {0, 0, 255}));
  connect(LP2.outlet, flowJoinTurbLP.in1) annotation(
    Line(points = {{76, 76}, {76, 36}, {34, 36}, {34, 26}}, color = {0, 0, 255}));
  connect(LP1.shaft_a, LP2.shaft_a) annotation(
    Line(points = {{13.2, 60}, {45.2, 60}}));
  connect(flowSplitTurbLP.out2, LP1.inlet) annotation(
    Line(points = {{26, 94}, {26, 76}, {16, 76}}, color = {0, 0, 255}));
  connect(LP1.outlet, flowJoinTurbLP.in2) annotation(
    Line(points = {{-16, 76}, {-16, 36}, {26, 36}, {26, 26}}, color = {0, 0, 255}));
  connect(HP.shaft_b, LP1.shaft_b) annotation(
    Line(points = {{-37.2, 60}, {-13.2, 60}}));
  connect(valveVapAdmissionTurb.outlet, HP.inlet) annotation(
    Line(points = {{-80, 76}, {-66, 76}}, color = {0, 0, 255}));
  connect(generator.port, powerSensorElectric.port_a) annotation(
    Line(points = {{168.6, 60}, {189.3, 60}, {189.3, -40}, {190, -40}}, color = {0, 0, 255}));
  connect(powerSensorElectric.port_b, load.port) annotation(
    Line(points = {{190, -60}, {190, -66.5}, {190.2, -66.5}, {190.2, -113}}, color = {0, 0, 255}));
  connect(powerSensorMech.flange_b, generator.shaft) annotation(
    Line(points = {{140, 60}, {152, 60}}));
  connect(HP.outlet, pressDropLin.inlet) annotation(
    Line(points = {{-34, 76}, {-34, 120}, {-20, 120}}, color = {0, 0, 255}));
  connect(pressDropLin.outlet, flowSplitTurbLP.in1) annotation(
    Line(points = {{0, 120}, {30, 120}, {30, 106}}, color = {0, 0, 255}));
  connect(steamGenBayonnet.outfl, SG1.infl) annotation(
    Line(points = {{-100, -100}, {-100, -200}, {-180, -200}, {-180, -180}}, color = {0, 0, 255}));
  connect(SG1.outfl, SG2.infl) annotation(
    Line(points = {{-180, -160}, {-180, -150}}, color = {0, 0, 255}));
  connect(SG2.outfl, SG3.infl) annotation(
    Line(points = {{-180, -130}, {-180, -120}}, color = {0, 0, 255}));
  connect(SG3.outfl, SG4.infl) annotation(
    Line(points = {{-180, -100}, {-180, -90}}, color = {0, 0, 255}));
  connect(SG4.outfl, SG5.infl) annotation(
    Line(points = {{-180, -70}, {-180, -60}}, color = {0, 0, 255}));
  connect(SG5.outfl, SG6.infl) annotation(
    Line(points = {{-180, -40}, {-180, -30}}, color = {0, 0, 255}));
  connect(SG6.outfl, SG7.infl) annotation(
    Line(points = {{-180, -10}, {-180, 0}}, color = {0, 0, 255}));
  connect(SG7.outfl, SG8.infl) annotation(
    Line(points = {{-180, 20}, {-180, 30}}, color = {0, 0, 255}));
  connect(heatSourceSG1.wall, SG1.wall) annotation(
    Line(points = {{-197, -170}, {-185, -170}}, color = {255, 127, 0}));
  connect(heatSourceSG2.wall, SG2.wall) annotation(
    Line(points = {{-197, -140}, {-185, -140}}, color = {255, 127, 0}));
  connect(heatSourceSG3.wall, SG3.wall) annotation(
    Line(points = {{-197, -110}, {-185, -110}}, color = {255, 127, 0}));
  connect(heatSourceSG4.wall, SG4.wall) annotation(
    Line(points = {{-197, -80}, {-185, -80}}, color = {255, 127, 0}));
  connect(heatSourceSG5.wall, SG5.wall) annotation(
    Line(points = {{-197, -50}, {-185, -50}}, color = {255, 127, 0}));
  connect(heatSourceSG6.wall, SG6.wall) annotation(
    Line(points = {{-197, -20}, {-185, -20}}, color = {255, 127, 0}));
  connect(heatSourceSG7.wall, SG7.wall) annotation(
    Line(points = {{-197, 10}, {-185, 10}}, color = {255, 127, 0}));
  connect(heatSourceSG8.wall, SG8.wall) annotation(
    Line(points = {{-197, 40}, {-185, 40}}, color = {255, 127, 0}));
  connect(steamGenPower1, signalProcess1.in_power) annotation(
    Line(points = {{-260, -170}, {-236, -170}}, color = {0, 0, 127}));
  connect(steamGenPower2, signalProcess2.in_power) annotation(
    Line(points = {{-260, -140}, {-236, -140}}, color = {0, 0, 127}));
  connect(steamGenPower3, signalProcess3.in_power) annotation(
    Line(points = {{-260, -110}, {-236, -110}}, color = {0, 0, 127}));
  connect(steamGenPower4, signalProcess4.in_power) annotation(
    Line(points = {{-260, -80}, {-236, -80}}, color = {0, 0, 127}));
  connect(steamGenPower5, signalProcess5.in_power) annotation(
    Line(points = {{-260, -50}, {-236, -50}}, color = {0, 0, 127}));
  connect(steamGenPower6, signalProcess6.in_power) annotation(
    Line(points = {{-260, -20}, {-236, -20}}, color = {0, 0, 127}));
  connect(steamGenPower7, signalProcess7.in_power) annotation(
    Line(points = {{-260, 10}, {-236, 10}}, color = {0, 0, 127}));
  connect(steamGenPower8, signalProcess8.in_power) annotation(
    Line(points = {{-260, 40}, {-236, 40}}, color = {0, 0, 127}));
  steamGenTemperature1 = SG1.T[2];
  steamGenTemperature2 = SG2.T[2];
  steamGenTemperature3 = SG3.T[2];
  steamGenTemperature4 = SG4.T[2];
  steamGenTemperature5 = SG5.T[2];
  steamGenTemperature6 = SG6.T[2];
  steamGenTemperature7 = SG7.T[2];
  steamGenTemperature8 = SG8.T[2];
  connect(signalProcess1.out_power, heatSourceSG1.power) annotation(
    Line(points = {{-215, -170}, {-205, -170}}, color = {0, 0, 127}));
  connect(signalProcess2.out_power, heatSourceSG2.power) annotation(
    Line(points = {{-215, -140}, {-205, -140}}, color = {0, 0, 127}));
  connect(signalProcess3.out_power, heatSourceSG3.power) annotation(
    Line(points = {{-215, -110}, {-205, -110}}, color = {0, 0, 127}));
  connect(signalProcess4.out_power, heatSourceSG4.power) annotation(
    Line(points = {{-215, -80}, {-205, -80}}, color = {0, 0, 127}));
  connect(signalProcess5.out_power, heatSourceSG5.power) annotation(
    Line(points = {{-215, -50}, {-205, -50}}, color = {0, 0, 127}));
  connect(signalProcess6.out_power, heatSourceSG6.power) annotation(
    Line(points = {{-215, -20}, {-205, -20}}, color = {0, 0, 127}));
  connect(signalProcess7.out_power, heatSourceSG7.power) annotation(
    Line(points = {{-215, 10}, {-205, 10}}, color = {0, 0, 127}));
  connect(signalProcess8.out_power, heatSourceSG8.power) annotation(
    Line(points = {{-215, 40}, {-205, 40}}, color = {0, 0, 127}));
  injectedPower = heatSourceSG1.power + heatSourceSG2.power + heatSourceSG3.power + heatSourceSG4.power + heatSourceSG5.power + heatSourceSG6.power + heatSourceSG7.power + heatSourceSG8.power;
  connect(flowJoin.in1, flowJoinTurbLP.out) annotation(
    Line(points = {{44, -6}, {30, -6}, {30, 14}}, color = {0, 0, 255}));
  connect(flowJoin.out, sinkPressure.flange) annotation(
    Line(points = {{56, -10}, {80, -10}}, color = {0, 0, 255}));
  connect(flowSplitBypass.in1, sensW.outlet) annotation(
    Line(points = {{-126, 72}, {-158, 72}}, color = {0, 0, 255}));
  connect(flowSplitBypass.out1, valveVapAdmissionTurb.inlet) annotation(
    Line(points = {{-114, 76}, {-100, 76}}, color = {0, 0, 255}));
  connect(valveVapBypass.outlet, flowJoin.in2) annotation(
    Line(points = {{-20, -14}, {44, -14}}, color = {0, 0, 255}));
  connect(flowSplitBypass.out2, valveVapBypass.inlet) annotation(
    Line(points = {{-114, 68}, {-100, 68}, {-100, -14}, {-40, -14}}, color = {0, 0, 255}));
  connect(valveLiqChoked.outlet, steamGenBayonnet.infl) annotation(
    Line(points = {{-80, -70}, {-100, -70}, {-100, -80}}, color = {0, 0, 255}));
  connect(frequencySensor.port, generator.port) annotation(
    Line(points = {{160, 150}, {174, 150}, {174, 60}, {168, 60}}, color = {0, 0, 255}));
  connect(inertia.flange_b, powerSensorMech.flange_a) annotation(
    Line(points = {{108, 60}, {120, 60}}));
  connect(LP2.shaft_b, inertia.flange_a) annotation(
    Line(points = {{72, 60}, {88, 60}}));
  connect(sensW.inlet, SG8.outfl) annotation(
    Line(points = {{-170, 72}, {-180, 72}, {-180, 50}}, color = {0, 0, 255}));
  connect(stepActivateExtReactivity.y, greaterThreshold.u) annotation(
    Line(points = {{261, -150}, {278, -150}}, color = {0, 0, 127}));
  connect(gain.u, powerSensorElectric.P) annotation(
    Line(points = {{248, -50}, {199, -50}}, color = {0, 0, 127}));
  connect(addPowers.y, corePowerMeasure.u) annotation(
    Line(points = {{271, 0}, {278, 0}}, color = {0, 0, 127}));
  connect(addPowers.u1, innerCorePower) annotation(
    Line(points = {{248, 6}, {239, 6}, {239, 20}, {220, 20}}, color = {0, 0, 127}));
  connect(addPowers.u2, outerCorePower) annotation(
    Line(points = {{248, -6}, {239, -6}, {239, -20}, {220, -20}}, color = {0, 0, 127}));
  connect(rampLoad.y, add3.u2) annotation(
    Line(points = {{131, -130}, {148, -130}}, color = {0, 0, 127}));
  connect(add3.y, load.referencePower) annotation(
    Line(points = {{171, -130}, {183, -130}}, color = {0, 0, 127}));
  connect(valveLiqChoked.inlet, sourcePressure.flange) annotation(
    Line(points = {{-60, -70}, {-28, -70}, {-28, -76}}, color = {0, 0, 255}));
  connect(rampLoad2.y, add3.u3) annotation(
    Line(points = {{131, -160}, {139, -160}, {139, -138}, {147, -138}}, color = {0, 0, 127}));
  connect(const.y, add3.u1) annotation(
    Line(points = {{131, -100}, {139, -100}, {139, -122}, {147, -122}}, color = {0, 0, 127}));
  connect(stepActivateExtSource.y, greaterThreshold1.u) annotation(
    Line(points = {{261, 90}, {278, 90}}, color = {0, 0, 127}));
  connect(gain.y, pidExternalSource.u_s) annotation(
    Line(points = {{272, -50}, {310, -50}, {310, 130}, {378, 130}}, color = {204, 0, 0}));
  connect(switch.u2, greaterThreshold1.y) annotation(
    Line(points = {{338, 90}, {301, 90}}, color = {255, 0, 255}));
  connect(gain.y, pidReactivity.u_s) annotation(
    Line(points = {{272, -50}, {310, -50}, {310, -110}, {378, -110}}, color = {204, 0, 0}));
  connect(greaterThreshold.y, switch1.u2) annotation(
    Line(points = {{301, -150}, {338, -150}}, color = {255, 0, 255}));
  connect(switch.y, pidExternalSource.u_m) annotation(
    Line(points = {{361, 90}, {390, 90}, {390, 118}}, color = {0, 0, 127}));
  connect(corePowerMeasure.y, switch.u3) annotation(
    Line(points = {{302, 0}, {330, 0}, {330, 82}, {338, 82}}, color = {0, 0, 127}));
  connect(gain.y, switch.u1) annotation(
    Line(points = {{272, -50}, {310, -50}, {310, 98}, {338, 98}}, color = {204, 0, 0}));
  connect(switch1.y, pidReactivity.u_m) annotation(
    Line(points = {{361, -150}, {390, -150}, {390, -122}}, color = {0, 0, 127}));
  connect(gain.y, switch1.u1) annotation(
    Line(points = {{272, -50}, {310, -50}, {310, -142}, {338, -142}}, color = {204, 0, 0}));
  connect(corePowerMeasure.y, switch1.u3) annotation(
    Line(points = {{302, 0}, {330, 0}, {330, -158}, {338, -158}}, color = {0, 0, 127}));
  connect(pidReactivity.y, externalReactivity) annotation(
    Line(points = {{402, -110}, {430, -110}}, color = {0, 0, 127}));
  connect(pidExternalSource.y, externalSourceModulation) annotation(
    Line(points = {{402, 130}, {430, 130}}, color = {0, 0, 127}));
  connect(frequencyOutput, frequencySensor.f) annotation(
    Line(points = {{120, 150}, {140, 150}}, color = {0, 0, 127}));
  connect(valveVapAdmissionTurb.theta, limiter.y) annotation(
    Line(points = {{-90, 84}, {-90, 180}, {39, 180}}, color = {0, 0, 127}));
  connect(add.y, limiter.u) annotation(
    Line(points = {{79, 180}, {62, 180}}, color = {0, 0, 127}));
  connect(add.u2, valveControlInput) annotation(
    Line(points = {{102, 174}, {102, 170}, {130, 170}}, color = {0, 0, 127}));
  connect(defaultValve.y, add.u1) annotation(
    Line(points = {{119, 190}, {102, 190}, {102, 186}}, color = {0, 0, 127}));
  annotation(
    uses(ThermoPower(version = "3.1"), Modelica(version = "3.2.3")),
    Diagram(coordinateSystem(extent = {{-280, 220}, {440, -200}}), graphics = {Rectangle(origin = {330, -140}, lineColor = {114, 159, 207}, fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-110, 50}, {110, -50}}), Rectangle(origin = {100, 180}, lineColor = {114, 159, 207}, fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-80, 40}, {80, -40}}), Text(origin = {395, -180}, extent = {{-45, 10}, {45, -10}}, textString = "External Reactivity 
Controller"), Text(origin = {50, 210}, extent = {{-30, 10}, {30, -10}}, textString = "Turbine-admission
Controller"), Rectangle(origin = {330, 100}, lineColor = {114, 159, 207}, fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-110, 50}, {110, -50}}), Text(origin = {395, 60}, extent = {{-45, 10}, {45, -10}}, textString = "External Source
Controller")}),
    Icon(coordinateSystem(extent = {{-200, -200}, {200, 200}}), graphics = {Rectangle(origin = {-160, -50}, fillColor = {114, 159, 207}, fillPattern = FillPattern.VerticalCylinder, extent = {{-40, 150}, {40, -150}}), Rectangle(fillColor = {186, 189, 182}, extent = {{-200, 200}, {200, -200}}), Text(origin = {0, -1}, extent = {{-80, 79}, {80, -79}}, textString = "Secondary
Circuit")}),
    version = "");
end SecondaryCircuitWithSimulink;