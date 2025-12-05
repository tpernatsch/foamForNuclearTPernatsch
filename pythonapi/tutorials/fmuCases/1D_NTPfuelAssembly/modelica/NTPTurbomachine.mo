model NTPTurbomachine
  extends Modelica.Icons.Example;
  package Medium = H2;
  parameter Real nozzleIR = 0.97;
  parameter Real nozzleOR = 1;
  parameter Real wnom = 30;
  ThermoPower.Gas.Compressor Compressor1(redeclare package Medium = Medium, pstart_in = 0.343e5, Tstart_in = 244.4, explicitIsentropicEnthalpy = true, Tstart_out = 600, pstart_out = 8.29e5, Ndesign = 523.3, Tdes_in = 244.4, Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix, tablePhic = tablePhicC, tableEta = tableEtaC, tablePR = tablePR) annotation(
    Placement(transformation(origin = {13.6, -22.4}, extent = {{-73.6, 38.4}, {-41.6, 6.4}})));
  ThermoPower.Gas.TurbineStodola Turbine1(redeclare package Medium = Medium, pstart_in = 7.85e5, pstart_out = 1.52e5, Tstart_out = 800, Tstart_in = 1390, Ndesign = 523.3, Tdes_in = 1400, fixedEta = false, wnom = wnom, Table = ThermoPower.Choices.TurboMachinery.TableTypes.matrix, tableEta = tableEtaT) annotation(
    Placement(transformation(origin = {-46.8, -22.4}, extent = {{92.8, 38.4}, {124.8, 6.4}})));
  ThermoPower.Gas.SourcePressure SourceP1(redeclare package Medium = Medium, T = 244.4, p0 = 0.3447e5) annotation(
    Placement(transformation(origin = {-60, 36}, extent = {{-100, -16}, {-80, 4}})));
  ThermoPower.Gas.SinkPressure SinkP1(redeclare package Medium = Medium, p0 = 1.52e5, T = 800, use_in_p0 = true, use_in_T = false) annotation(
    Placement(transformation(origin = {38, -6}, extent = {{82, -16}, {102, 4}})));
  Modelica.Mechanics.Rotational.Components.Inertia Inertia1(J = 50) annotation(
    Placement(transformation(origin = {-6, 14}, extent = {{6, -24}, {26, -4}})));
  ThermoPower.Gas.PressDrop PressDrop3(redeclare package Medium = Medium, FFtype = ThermoPower.Choices.PressDrop.FFtypes.OpPoint, A = 1, wnom = wnom, rhonom = 2, dpnom = 26000, pstart = 811000, Tstart = 1370) annotation(
    Placement(transformation(origin = {0, -144}, extent = {{34, 0}, {54, 20}}, rotation = 90)));
  ThermoPower.Gas.PressDrop PressDrop2(FFtype = ThermoPower.Choices.PressDrop.FFtypes.OpPoint, A = 1, redeclare package Medium = Medium, wnom = wnom, rhonom = 4.7, dpnom = 18000, pstart = 829000, Tstart = 600) annotation(
    Placement(transformation(origin = {-40, -58}, extent = {{-20, 0}, {0, 20}}, rotation = -90)));
  ThermoPower.Gas.PressDrop PressDrop1(FFtype = ThermoPower.Choices.PressDrop.FFtypes.OpPoint, A = 1, redeclare package Medium = Medium, wnom = wnom, rhonom = 0.48, dpnom = 170, pstart = 34470, Tstart = 244.4) annotation(
    Placement(transformation(origin = {-64, -52}, extent = {{-72, -16}, {-52, 4}}, rotation = -90)));
  inner ThermoPower.System system annotation(
    Placement(transformation(origin = {200, -340}, extent = {{80, 80}, {100, 100}})));
  parameter Real Dhyd = 0.01;
  ThermoPower.Gas.Flow1DFV NozzleWall(redeclare package Medium = Medium, L = 2, Dhyd = 4*3.14159265*(nozzleOR*nozzleOR - nozzleIR*nozzleIR)/(2*3.141592*nozzleOR), A = 3.14159265*(nozzleOR*nozzleOR - nozzleIR*nozzleIR), omega = 2*3.14159265*1, wnom = wnom, N = 11, pstart = 811000, Tstartbar = 873.15, Nt = 1, Tstartin = 600, Tstartout = 1370, QuasiStatic = true) annotation(
    Placement(transformation(origin = {-10, -150}, extent = {{-10, 10}, {10, -10}}, rotation = 90)));
  /*ThermoPower.Thermal.TempSource1DFV tempSource1DFV(Nw = 10) annotation(
                                                                                                    Placement(transformation(origin = {20, 50}, extent = {{-10, -10}, {10, 10}})));
                                                                                                  Modelica.Blocks.Sources.Constant const(k = 350) annotation(
                                                                                                    Placement(transformation(origin = {-10, 90}, extent = {{-10, -10}, {10, 10}})));*/
  ThermoPower.Gas.SinkPressure sinkPressure(redeclare package Medium = Medium, p0 = 30e5/*4.71*/, T = 573.15) annotation(
    Placement(transformation(origin = {70, -230}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  ThermoPower.Gas.Flow1DFV NozzleChamber(A = 3.14159265*nozzleIR*nozzleIR, Dhyd = 2*nozzleIR, L = 2, redeclare package Medium = Medium, N = 11, Nt = 1, QuasiStatic = true, Tstartbar(displayUnit = "K") = 300, omega = 2*3.14159265*nozzleIR, pstart = 3e6, wnom = wnom) annotation(
    Placement(transformation(origin = {70, -150}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
  ThermoPower.Thermal.MetalTubeFV metalTubeFV1(rint = nozzleIR, rext = nozzleOR, Nw = 10, L = 2, rhomcm = 0.49e3*7680, lambda = 45) annotation(
    Placement(transformation(origin = {10, -150}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  ThermoPower.Thermal.CounterCurrentFV counterCurrentFV(Nw = 10) annotation(
    Placement(transformation(origin = {30, -150}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  ThermoPower.Gas.SourceMassFlow sourceMassFlow(p0 = 541000, T(displayUnit = "K") = 300, redeclare package Medium = Medium, w0 = wnom, use_in_T = true, use_in_w0 = true) annotation(
    Placement(transformation(origin = {70, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  ThermoPower.Gas.SensT1 sensT2(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {94, 74}, extent = {{-10, -10}, {10, 10}})));
  ThermoPower.Gas.SensW sensW(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {100, -8}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = 0.5, y(start = 526.85), y_start = 526.85) annotation(
    Placement(transformation(origin = {190, 80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder1(T = 1, k = 1/1500) annotation(
    Placement(transformation(origin = {190, 40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput Tturbine_out annotation(
    Placement(transformation(origin = {230, 80}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {210, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput mFlowTurbine_out annotation(
    Placement(transformation(origin = {230, 40}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {210, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput pTurbine_in annotation(
    Placement(transformation(origin = {240, 10}, extent = {{20, -20}, {-20, 20}}), iconTransformation(origin = {-220, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder2(T = 1) annotation(
    Placement(transformation(origin = {200, 10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = 45e5, uMin = 1.5e5) annotation(
    Placement(transformation(origin = {170, 10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Interfaces.RealInput Tnozzle_in annotation(
    Placement(transformation(origin = {280, -110}, extent = {{20, -20}, {-20, 20}}), iconTransformation(origin = {-220, -82}, extent = {{-20, -20}, {20, 20}})));
  ThermoPower.Thermal.MetalTubeFV metalTubeFV11(L = 2, Nw = 10, lambda = 45, rext = nozzleOR, rhomcm = 0.49e3*7680, rint = nozzleIR) annotation(
    Placement(transformation(origin = {50, -150}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder3(T = 1) annotation(
    Placement(transformation(origin = {190, -110}, extent = {{10, -10}, {-10, 10}})));
  ThermoPower.Gas.SensP sensP(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {110, -144}, extent = {{-10, 10}, {10, -10}})));
  ThermoPower.Gas.PressDrop pressDrop(wnom = wnom, dpnom = 1e4, Kf = 0.1, redeclare package Medium = Medium, FFtype = ThermoPower.Choices.PressDrop.FFtypes.Kinetic, K = 10, A = 1, pstart = 3e6) annotation(
    Placement(transformation(origin = {70, -190}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput pNozzleChamber_out annotation(
    Placement(transformation(origin = {230, -150}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {210, -70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder4(T = 1) annotation(
    Placement(transformation(origin = {190, -150}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Nonlinear.Limiter limiter2(uMax = 60e5, uMin = 20e5) annotation(
    Placement(transformation(origin = {150, -150}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Nonlinear.Limiter limiter1(uMax = 3000, uMin = 300) annotation(
    Placement(transformation(origin = {150, -110}, extent = {{10, -10}, {-10, 10}})));
  ThermoPower.Gas.SensT1 sensTnozzleWall_out(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {10, -120}, extent = {{-10, -10}, {10, 10}})));
  ThermoPower.Gas.SensT1 sensTnozzleWall_in(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {10, -190}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput corePower_in annotation(
    Placement(transformation(origin = {-220, -40}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-220, -160}, extent = {{-20, -20}, {20, 20}})));
  ControlSystem controlSystem(wnom = wnom)  annotation(
    Placement(transformation(origin = {-160, -40}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealInput mFlowCoreOutlet_in annotation(
    Placement(transformation(origin = {240, -70}, extent = {{20, -20}, {-20, 20}}), iconTransformation(origin = {-220, 80}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder5(T = 1) annotation(
    Placement(transformation(origin = {190, -70}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Nonlinear.Limiter limiter3(uMax = 2*wnom, uMin = 0.2*wnom) annotation(
    Placement(transformation(origin = {150, -70}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  ThermoPower.Gas.ValveLin valveLin(Kv = 0.01, redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {-90, 30}, extent = {{-10, 10}, {10, -10}}, rotation = -0)));
  Modelica.Blocks.Continuous.FirstOrder firstOrder6(T = 1, initType = Modelica.Blocks.Types.Init.InitialOutput, y_start = 0.5) annotation(
    Placement(transformation(origin = {-110, -30}, extent = {{-10, -10}, {10, 10}})));
protected
  /*
  parameter Real tableEtaC[6, 4] = [
    0, 95,      100,      105;
    1, 82.5e-2, 81e-2,    80.5e-2;
    2, 84e-2,   82.9e-2,  82e-2;
    3, 83.2e-2, 82.2e-2,  81.5e-2;
    4, 82.5e-2, 81.2e-2,  79e-2;
    5, 79.5e-2, 78e-2,    76.5e-2
  ];
  parameter Real tablePhicC[6, 4] = [
    0, 95,          100,         105;
    1, 38.3e-3*0.3, 43e-3*0.3,   46.8e-3*0.3;
    2, 39.3e-3*0.3, 43.8e-3*0.3, 47.9e-3*0.3;
    3, 40.6e-3*0.3, 45.2e-3*0.3, 48.4e-3*0.3;
    4, 41.6e-3*0.3, 46.1e-3*0.3, 48.9e-3*0.3;
    5, 42.3e-3*0.3, 46.6e-3*0.3, 49.3e-3*0.3
  ];
  parameter Real tablePR[6, 4] = [
    0, 95, 100, 105;
    1, 22.6, 27, 32;
    2, 22, 26.6, 30.8;
    3, 20.8, 25.5, 29;
    4, 19, 24.3, 27.1;
    5, 17, 21.5, 24.2
  ];
  parameter Real tableEtaT[5, 4] = [
    1, 90, 100, 110;
    2.36, 89e-2, 89.5e-2, 89.3e-2;
    2.88, 90e-2, 90.6e-2, 90.5e-2;
    3.56, 90.5e-2, 90.6e-2, 90.5e-2;
    4.46, 90.2e-2, 90.3e-2, 90e-2
  ];
  // */
  // /*
  parameter Real tableEtaC[6, 5] = [
    0, 50,   95,      100,      105;
    1, 0.70, 82.5e-2, 81e-2,    80.5e-2;
    2, 0.73, 84e-2,   82.9e-2,  82e-2;
    3, 0.72, 83.2e-2, 82.2e-2,  81.5e-2;
    4, 0.70, 82.5e-2, 81.2e-2,  79e-2;
    5, 0.66, 79.5e-2, 78e-2,    76.5e-2
  ];
  parameter Real tablePhicC[6, 5] = [
    0, 50,        95,          100,         105;
    1, 0.020*0.3, 38.3e-3*0.3, 43e-3*0.3,   46.8e-3*0.3;
    2, 0.021*0.3, 39.3e-3*0.3, 43.8e-3*0.3, 47.9e-3*0.3;
    3, 0.022*0.3, 40.6e-3*0.3, 45.2e-3*0.3, 48.4e-3*0.3;
    4, 0.023*0.3, 41.6e-3*0.3, 46.1e-3*0.3, 48.9e-3*0.3;
    5, 0.024*0.3, 42.3e-3*0.3, 46.6e-3*0.3, 49.3e-3*0.3
  ];
  parameter Real tablePR[6, 5] = [
    0, 50,  95,    100,  105;
    1, 5.0, 22.6,  27,   32;
    2, 6.0, 22,    26.6, 30.8;
    3, 7.5, 20.8,  25.5, 29;
    4, 8.5, 19,    24.3, 27.1;
    5, 9.0, 17,    21.5, 24.2
  ];
  parameter Real tableEtaT[5, 17] = [
    1,     1.0,   10,    20,    30,    40,    50,    60,    70,    80,    90,     100,    110,    120,    130,    140,    150;
    2.36,  0.30,  0.45,  0.60,  0.72,  0.80,  0.85,  0.87,  0.88,  0.885, 0.89,   0.895,  0.893,  0.885,  0.87,   0.85,   0.82;
    2.88,  0.35,  0.50,  0.65,  0.77,  0.84,  0.88,  0.89,  0.895, 0.898, 0.90,   0.906,  0.905,  0.895,  0.88,   0.86,   0.83;
    3.56,  0.38,  0.52,  0.67,  0.79,  0.86,  0.90,  0.905, 0.906, 0.905, 0.905,  0.906,  0.905,  0.895,  0.88,   0.86,   0.83;
    4.46,  0.36,  0.50,  0.65,  0.76,  0.83,  0.88,  0.895, 0.898, 0.90,  0.902,  0.903,  0.900,  0.89,   0.87,   0.84,   0.81
  ];
  // */
equation
  connect(Compressor1.shaft_b, Inertia1.flange_a) annotation(
    Line(points = {{-34, 0}, {0, 0}}, thickness = 0.5));
  connect(Inertia1.flange_b, Turbine1.shaft_a) annotation(
    Line(points = {{20, 0}, {52, 0}}, thickness = 0.5));
  connect(Compressor1.outlet, PressDrop2.inlet) annotation(
    Line(points = {{-31, -13}, {-30, -13}, {-30, -38}}, color = {159, 159, 223}, thickness = 0.5));
  connect(PressDrop1.outlet, Compressor1.inlet) annotation(
    Line(points = {{-70, 0}, {-70, -13}, {-57, -13}}, color = {159, 159, 223}, thickness = 0.5));
  connect(PressDrop2.outlet, NozzleWall.infl) annotation(
    Line(points = {{-30, -58}, {-30, -172}, {-10, -172}, {-10, -160}}, color = {159, 159, 223}, thickness = 0.5));
  connect(NozzleWall.outfl, PressDrop3.inlet) annotation(
    Line(points = {{-10, -140}, {-10, -110}}, color = {159, 159, 223}, thickness = 0.5));
  connect(metalTubeFV1.int, counterCurrentFV.side2) annotation(
    Line(points = {{13, -150}, {27, -150}}, color = {255, 127, 0}));
  connect(metalTubeFV1.ext, NozzleWall.wall) annotation(
    Line(points = {{6.9, -150}, {-5.1, -150}}, color = {255, 127, 0}));
  connect(metalTubeFV11.int, NozzleChamber.wall) annotation(
    Line(points = {{53, -150}, {65, -150}}, color = {255, 127, 0}));
  connect(metalTubeFV11.ext, counterCurrentFV.side1) annotation(
    Line(points = {{46.9, -150}, {32.9, -150}}, color = {255, 127, 0}));
  connect(sourceMassFlow.flange, NozzleChamber.infl) annotation(
    Line(points = {{70, -120}, {70, -140}}, color = {159, 159, 223}, thickness = 0.5));
  connect(firstOrder.u, sensT2.T) annotation(
    Line(points = {{178, 80}, {101, 80}}, color = {0, 0, 127}));
  connect(firstOrder1.u, sensW.w) annotation(
    Line(points = {{178, 40}, {115, 40}, {115, -2}, {107, -2}}, color = {0, 0, 127}));
  connect(firstOrder.y, Tturbine_out) annotation(
    Line(points = {{201, 80}, {230, 80}}, color = {0, 0, 127}));
  connect(firstOrder1.y, mFlowTurbine_out) annotation(
    Line(points = {{201, 40}, {230, 40}}, color = {0, 0, 127}));
  connect(pTurbine_in, firstOrder2.u) annotation(
    Line(points = {{240, 10}, {212, 10}}, color = {0, 0, 127}));
  connect(firstOrder2.y, limiter.u) annotation(
    Line(points = {{189, 10}, {182, 10}}, color = {0, 0, 127}));
  connect(limiter.y, SinkP1.in_p0) annotation(
    Line(points = {{159, 10}, {124, 10}, {124, -6}}, color = {0, 0, 127}));
  connect(sensW.outlet, SinkP1.flange) annotation(
    Line(points = {{106, -12}, {120, -12}}, color = {159, 159, 223}, thickness = 0.5));
  connect(Tnozzle_in, firstOrder3.u) annotation(
    Line(points = {{280, -110}, {202, -110}}, color = {0, 0, 127}, thickness = 0.5));
  connect(sensP.flange, sourceMassFlow.flange) annotation(
    Line(points = {{110, -140}, {110, -130}, {70, -130}, {70, -120}}, color = {159, 159, 223}, thickness = 0.5));
  connect(NozzleChamber.outfl, pressDrop.inlet) annotation(
    Line(points = {{70, -160}, {70, -180}}, color = {159, 159, 223}, thickness = 0.5));
  connect(pressDrop.outlet, sinkPressure.flange) annotation(
    Line(points = {{70, -200}, {70, -220}}, color = {159, 159, 223}, thickness = 0.5));
  connect(sourceMassFlow.in_T, limiter1.y) annotation(
    Line(points = {{75, -110}, {139, -110}}, color = {0, 0, 127}));
  connect(limiter1.u, firstOrder3.y) annotation(
    Line(points = {{162, -110}, {179, -110}}, color = {0, 0, 127}));
  connect(sensP.p, limiter2.u) annotation(
    Line(points = {{117, -150}, {137, -150}}, color = {0, 0, 127}));
  connect(limiter2.y, firstOrder4.u) annotation(
    Line(points = {{161, -150}, {177, -150}}, color = {0, 0, 127}));
  connect(firstOrder4.y, pNozzleChamber_out) annotation(
    Line(points = {{201, -150}, {229, -150}}, color = {0, 0, 127}));
  connect(NozzleWall.outfl, sensTnozzleWall_out.flange) annotation(
    Line(points = {{-10, -140}, {-10, -124}, {10, -124}}, color = {159, 159, 223}, thickness = 0.5));
  connect(NozzleWall.infl, sensTnozzleWall_in.flange) annotation(
    Line(points = {{-10, -160}, {-10, -194}, {10, -194}}, color = {159, 159, 223}, thickness = 0.5));
  connect(corePower_in, controlSystem.power) annotation(
    Line(points = {{-220, -40}, {-182, -40}}, color = {0, 0, 127}));
  connect(limiter3.y, sourceMassFlow.in_w0) annotation(
    Line(points = {{140, -70}, {100, -70}, {100, -104}, {76, -104}}, color = {0, 0, 127}));
  connect(limiter3.u, firstOrder5.y) annotation(
    Line(points = {{162, -70}, {180, -70}}, color = {0, 0, 127}));
  connect(firstOrder5.u, mFlowCoreOutlet_in) annotation(
    Line(points = {{202, -70}, {240, -70}}, color = {0, 0, 127}));
  connect(sensT2.flange, sensW.inlet) annotation(
    Line(points = {{94, 70}, {94, -12}}, color = {159, 159, 223}, thickness = 0.5));
  connect(Turbine1.outlet, sensW.inlet) annotation(
    Line(points = {{75, -13}, {75, -12}, {94, -12}}, color = {159, 159, 223}, thickness = 0.5));
  connect(controlSystem.turbineValveOpening, firstOrder6.u) annotation(
    Line(points = {{-139, -30}, {-123, -30}}, color = {0, 0, 127}));
  connect(sensW.w, controlSystem.massFlow) annotation(
    Line(points = {{107, -2}, {115, -2}, {115, 52}, {-200, 52}, {-200, -30}, {-182, -30}}, color = {0, 0, 127}));
  connect(firstOrder6.y, valveLin.cmd) annotation(
    Line(points = {{-98, -30}, {-90, -30}, {-90, 23}}, color = {0, 0, 127}));
  connect(Tnozzle_in, controlSystem.Tnozzle) annotation(
    Line(points = {{280, -110}, {240, -110}, {240, -254}, {-200, -254}, {-200, -50}, {-182, -50}}, color = {0, 0, 127}, thickness = 0.5));
  connect(PressDrop3.outlet, Turbine1.inlet) annotation(
    Line(points = {{-10, -90}, {-10, -50}, {50, -50}, {50, -12}}, color = {159, 159, 223}));
  connect(SourceP1.flange, valveLin.inlet) annotation(
    Line(points = {{-140, 30}, {-100, 30}}, color = {159, 159, 223}));
  connect(valveLin.outlet, PressDrop1.inlet) annotation(
    Line(points = {{-80, 30}, {-70, 30}, {-70, 20}}, color = {159, 159, 223}));
initial equation
  Inertia1.phi = 0;
  Inertia1.w = 523 + 10.2598 - 1.22443e-5;
  der(Inertia1.w) = 0;
  annotation(
    Documentation(info = "<html>
This is a simplified turbomachine model for NTP.
</html>"),
    experiment(StopTime = 1000, Tolerance = 1e-06, StartTime = 0, Interval = 2),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-240, 100}, {300, -260}}, initialScale = 0.1), graphics = {Rectangle(origin = {35, -164}, lineColor = {53, 132, 228}, lineThickness = 1, extent = {{-57, 84}, {57, -84}}), Text(origin = {-10, -240}, extent = {{-10, 6}, {10, -6}}, textString = "Nozzle")}),
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -200}, {200, 200}}, initialScale = 0.1), graphics = {Rectangle(lineColor = {170, 170, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-200, 200}, {200, -200}}), Text(textColor = {170, 170, 255}, extent = {{-140, 140}, {140, -140}}, textString = "P")}),
    version = "",
    uses(ThermoPower(version = "3.1"), Modelica(version = "3.2.3")),
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=uode --allowNonStandardModelica=reinitInAlgorithms -d=initialization,NLSanalyticJacobian");
end NTPTurbomachine;