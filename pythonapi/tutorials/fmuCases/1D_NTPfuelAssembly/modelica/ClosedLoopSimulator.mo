model ClosedLoopSimulator
  extends Modelica.Icons.Example;
  /*
    ThermoPower.Examples.RankineCycle.Models.PID
                     pID(
      Ti=5,
      PVmin=0,
      PVmax=2,
      CSmin=0,
      CSmax=4,
      steadyStateInit=true,
      Kp= 1,
      holdWhenSimplified=true, PVstart = 1, CSstart = 1) annotation (Placement(transformation(extent = {{-32, -10}, {-12, 10}})));
    */
  inner ThermoPower.System system annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  Modelica.Blocks.Continuous.LimPID PID(controllerType = Modelica.Blocks.Types.SimpleController.PI, Ti = 0.01, yMax = 10e5, k = 1, initType = Modelica.Blocks.Types.InitPID.InitialOutput, y_start = 1.52e5, yMin = 0)  annotation(
    Placement(transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}})));
  NTPTurbomachine2 nTPTurbomachine2 annotation(
    Placement(transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(duration = 20, height = 100, offset = 300, startTime = 200) annotation(
    Placement(transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(ramp.y, PID.u_s) annotation(
    Line(points = {{-58, 0}, {-42, 0}}, color = {0, 0, 127}));
  connect(PID.u_m, nTPTurbomachine2.Tnozzle) annotation(
    Line(points = {{-30, -12}, {-30, -40}, {60, -40}, {60, 0}, {40, 0}}, color = {0, 0, 127}));
  connect(PID.y, nTPTurbomachine2.p_out) annotation(
    Line(points = {{-18, 0}, {20, 0}}, color = {0, 0, 127}));
  annotation(
    Diagram(graphics),
    experiment(StopTime = 1000, Tolerance = 1e-006),
  __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=uode --allowNonStandardModelica=reinitInAlgorithms -d=initialization,NLSanalyticJacobian",
    Documentation(revisions = "<html>
<ul>
<li><i>10 Dec 2008</i>
by <a>Luca Savoldelli</a>:<br>
   First release.</li>
</ul>
</html>", info = "<html>
<p>This model simulates a simple continuous-time control system for the steam power plant. The generated power is controlled to the set point by a PI controller with anti-windup.</p>
<p>The model starts at steady state.
</html>"),
    uses(Modelica(version = "3.2.3")));
end ClosedLoopSimulator;