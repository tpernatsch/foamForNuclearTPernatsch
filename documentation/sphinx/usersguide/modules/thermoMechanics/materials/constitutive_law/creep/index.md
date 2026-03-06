# <b>Creep Models</b>

Creep models provide the creep strain rate (or equivalent creep contribution)
used by constitutive laws that include time-dependent inelastic deformation.

Creep models are selected where required by the chosen constitutive law. The
specific keyword and expected inputs are documented in the resepective creep law
page.

## <b>Classes</b>

Currently, OFFBEAT supports the following creep models:

- [fromLatestTime](../../../../classes/rheology/constitutiveLaws/creepModels/creepModel.md)  
  base creep model (also used to read creep-related state from the latest time).

- [powerLaw](../../../../classes/rheology/constitutiveLaws/creepModels/powerLaw.md)  
  power-law creep model.

- [constantPrincipalStress](../../../../classes/rheology/constitutiveLaws/creepModels/constantCreepPrincipalStress.md)  
  creep model based on a constant principal stress assumption.

- [correlationPrincipalStress](../../../../classes/rheology/constitutiveLaws/creepModels/correlationCreepPrincipalStress.md)  
  creep model based on a principal stress correlation.

- [UO2Matpro](../../../../classes/rheology/constitutiveLaws/creepModels/MatproCreepModel.md)  
  MATPRO-based creep model for UO₂ fuel.

- [UPuO2Malygin](../../../../classes/rheology/constitutiveLaws/creepModels/MalyginMOXCreepModel.md)  
  Malygin creep model for MOX fuel.

- [UPuO2Routbort](../../../../classes/rheology/constitutiveLaws/creepModels/RoutbortFastMOXCreepModel.md)  
  Routbort creep model for fast-reactor MOX fuel.

- [ZircaloyLimback](../../../../classes/rheology/constitutiveLaws/creepModels/ZircaloyLimback.md)  
  Limback creep model for Zircaloy.

- [ZircaloyLimbackLoca](../../../../classes/rheology/constitutiveLaws/creepModels/ZircaloyLimbackLoca.md)  
  Limback LOCA creep model for Zircaloy.

- [Steel1515TiTobbe](../../../../classes/rheology/constitutiveLaws/creepModels/TobbeCreepModel.md)  
  Tobbe creep model for 15-15Ti steel.

- [Steel1515TiAim1Tobbe](../../../../classes/rheology/constitutiveLaws/creepModels/Steel1515TiAim1Tobbe.md)  
  AIM1 + Tobbe creep model for 15-15Ti steel.

- [Steel1515TiDINTobbe](../../../../classes/rheology/constitutiveLaws/creepModels/TobbeDINCreepModel.md)  
  DIN + Tobbe creep model for 15-15Ti steel.

- [HastelloyZhang](../../../../classes/rheology/constitutiveLaws/creepModels/ZhangHastelloyCreepModel.md)  
  Zhang creep model for Hastelloy.

- [SiCMonolithic](../../../../classes/rheology/constitutiveLaws/creepModels/MonolithicSiCCreepModel.md)  
  creep model for monolithic SiC.

- [BufferParfume](../../../../classes/rheology/constitutiveLaws/creepModels/ParfumeBufferCreepModel.md)  
  PARFUME-based creep model for buffer layers.

- [PyCParfume](../../../../classes/rheology/constitutiveLaws/creepModels/ParfumePyCCreepModel.md)  
  PARFUME-based creep model for PyC layers.
