Motivation and Context
=====================

The evolution of nuclear fuel during irradiation is governed by complex
phenomena spanning multiple spatial and temporal scales. Most traditional
nuclear fuel performance codes have addressed this complexity by relying on
reduced-dimensional representations, such as one-dimensional or
quasi-two-dimensional models.

These approaches are well suited for many engineering-scale applications,
but they inherently limit the analysis of local effects, non-axisymmetric
conditions, and strongly spatial phenomena. In addition, many tools have
historically focused either on steady-state base irradiation or on specific
transient and accident scenarios, rather than supporting both within a
single modelling approach.

Over the first two decades of the 2000s, the fuel performance community has
shown an increasing interest in higher-fidelity analysis tools capable of
addressing multi-dimensional, multi-physics, and multi-scale problems. This
trend reflects the need to better understand local mechanisms affecting fuel
integrity and to support the interpretation of experimental observations
that cannot be explained using reduced-order models alone.

The development of OFFBEAT was initiated in this context. A specific
motivation was provided by a fuel failure that occurred in a Swiss boiling
water reactor, which triggered an extended root cause analysis at PSI
:cite:`CliffordLocalPeakingBWR` aimed at identifying potential local and
asymmetric effects. This effort highlighted the need for a modelling tool
capable of resolving complex thermo-mechanical behaviour beyond the
assumptions typically adopted in traditional fuel performance analyses.

OFFBEAT was therefore conceived as a fully multi-dimensional fuel
performance code, with the objective of improving the understanding of
local fuel behaviour mechanisms and assessing their impact on fuel
integrity. The initial development of the code was driven primarily by the
work of Scolaro :cite:`ScolaroThesis`, :cite:`ScolaroOffbeatNED`, building on
the works of Jasak and Weller :cite:`JasakWellerLinearElasticity`, Tuković
:cite:`TukovicFluidSolid`, Cardiff :cite:`Cardiff30`, and Clifford
:cite:`CliffordLocalPeakingBWR`. Subsequent extensions and validation
activities have been carried out by multiple contributors, further
expanding the scope and capabilities of the code.

While early applications focused on the analysis of asymmetric heat
transfer conditions, the scope of OFFBEAT has since expanded to support a
broad range of steady-state and transient fuel behaviour studies.
