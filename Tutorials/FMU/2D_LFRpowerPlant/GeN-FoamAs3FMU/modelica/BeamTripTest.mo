model BeamTripTest
  Modelica.Blocks.Noise.UniformNoise uniformNoise(samplePeriod = 1, y_max = 1, y_min = 0)  annotation(
    Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = 1 - 7/(24*3600))  annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(uniformNoise.y, greaterThreshold.u) annotation(
    Line(points = {{-49, 0}, {-13, 0}}, color = {0, 0, 127}));

annotation(
    uses(Modelica(version = "3.2.3")),
  Documentation(info = "<html><head></head><body>Simple model of 1s beam trips that occurs 7 times per day in average.</body></html>"));
end BeamTripTest;