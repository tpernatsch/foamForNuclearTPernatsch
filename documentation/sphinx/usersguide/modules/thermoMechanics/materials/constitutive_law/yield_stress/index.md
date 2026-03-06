# <b>Yield Stress Models</b>

Yield stress models provide the yield strength used by plasticity-based
constitutive laws. They typically define the dependence of yield stress on
temperature, irradiation, strain, strain rate, or other state variables,
depending on the model.

Yield stress models are selected where required by the chosen constitutive law.
The specific keyword and expected inputs are documented in the constitutive law
page.

## <b>Classes</b>

Currently, OFFBEAT supports the following yield stress models:

- [constant](../../../../classes/rheology/constitutiveLaws/yieldStressModels/yieldStressConstant.md)  
  assigns a constant yield stress.

- [Fraptran](../../../../classes/rheology/constitutiveLaws/yieldStressModels/yieldStressFraptran.md)  
  FRAPTRAN-based yield stress model.

- [hardening](../../../../classes/rheology/constitutiveLaws/yieldStressModels/hardening.md)  
  yield stress model including hardening effects.
