model ExternalReactivityController
  Modelica.Blocks.Sources.Ramp rampReactivity(duration = 1, height = 0.0005767131, offset = 0, startTime = 100)  annotation(
    Placement(visible = true, transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput externalReactivity annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(rampReactivity.y, externalReactivity);
  annotation(
    uses(Modelica(version = "3.2.3")));
end ExternalReactivityController;