model OpenLoopSimulator_ControlSystem
  extends Modelica.Icons.Example;
  inner ThermoPower.System system annotation(
    Placement(transformation(extent = {{80, 80}, {100, 100}})));
  ControlSystem controlSystem annotation(
    Placement(transformation(extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Sources.Ramp ramp(height = 1e6, duration = 100, offset = 1e6, startTime = 100) annotation(
    Placement(transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 100)  annotation(
    Placement(transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = 10)  annotation(
    Placement(transformation(origin = {-56, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(const.y, controlSystem.massFlow) annotation(
    Line(points = {{-78, 50}, {-40, 50}, {-40, 10}, {-22, 10}}, color = {0, 0, 127}));
  connect(ramp.y, firstOrder.u) annotation(
    Line(points = {{-78, 0}, {-68, 0}}, color = {0, 0, 127}));
  connect(firstOrder.y, controlSystem.power) annotation(
    Line(points = {{-44, 0}, {-22, 0}}, color = {0, 0, 127}));
  annotation(
    Diagram(graphics),
    experiment(StopTime = 1000, Tolerance = 1e-06, StartTime = 0, Interval = 0.2),
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=uode --allowNonStandardModelica=reinitInAlgorithms -d=initialization,NLSanalyticJacobian",
    Documentation(revisions = "<html>
<ul>
<li><i>10 Dec 2008</i>
by <a>Luca Savoldelli</a>:<br>
   First release.</li>
</ul>
</html>", info = "<html>
<p>This model allows to simulate an open loop transients.
</html>"),
    __Dymola_experimentSetupOutput,
    uses(Modelica(version = "3.2.3")));
end OpenLoopSimulator_ControlSystem;