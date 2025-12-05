package H2 "Hydrogen model"
  extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
    mediumName="H2",
    data={Modelica.Media.IdealGases.Common.SingleGasesData.H2},
    fluidConstants={
          Modelica.Media.IdealGases.Common.FluidData.H2},
    substanceNames={"Hydrogen"},
    reference_X={1},
    referenceChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt25C);
end H2;
