model OpenLoopSimulator
  extends Modelica.Icons.Example;
  inner ThermoPower.System system
    annotation (Placement(transformation(extent={{80,80},{100,100}})));
  Modelica.Blocks.Sources.Ramp rampP(duration = 100, height = 30e5, offset = 1.52e5, startTime = 200) annotation(
    Placement(transformation(origin = {-30, -10}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp rampT(duration = 100, height = 2000, offset = 300, startTime = 500) annotation(
    Placement(transformation(origin = {-30, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant mFlow(k = 100)  annotation(
    Placement(transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant power(k = 900e6/1500)  annotation(
    Placement(transformation(origin = {-58, -92}, extent = {{-10, -10}, {10, 10}})));
  NTPTurbomachine3 nTPTurbomachine3 annotation(
    Placement(transformation(origin = {40, 0}, extent = {{-20, -20}, {20, 20}})));
equation
  connect(mFlow.y, nTPTurbomachine3.mFlowCoreOutlet_in) annotation(
    Line(points = {{-38, 50}, {0, 50}, {0, 8}, {18, 8}}, color = {0, 0, 127}));
  connect(rampP.y, nTPTurbomachine3.pTurbine_in) annotation(
    Line(points = {{-18, -10}, {-6, -10}, {-6, 0}, {18, 0}}, color = {0, 0, 127}));
  connect(rampT.y, nTPTurbomachine3.Tnozzle_in) annotation(
    Line(points = {{-18, -50}, {-2, -50}, {-2, -8}, {18, -8}}, color = {0, 0, 127}));
  connect(power.y, nTPTurbomachine3.corePower_in) annotation(
    Line(points = {{-46, -92}, {4, -92}, {4, -16}, {18, -16}}, color = {0, 0, 127}));
  annotation(
    Diagram(graphics),
    experiment(StopTime = 1000, __Dymola_NumberOfIntervals = 5000, Tolerance = 1e-006),
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
end OpenLoopSimulator;