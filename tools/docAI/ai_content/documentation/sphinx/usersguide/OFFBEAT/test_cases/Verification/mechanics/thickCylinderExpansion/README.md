# Expanding Cylinder

An infinite hollow cylinder initially has inner and outer radii of 1 cm and 2 cm, respectively. The inner surface expands radially outward over 75 s at a rate of approximately 1 mm/s, reaching a final inner radius of approximately 8.5 cm. The bulk and shear moduli of the material are 3.8 GPa and 40 GPa, respectively. The yield stress is \( \sigma_y = 500 \) MPa. We assume perfect plasticity.

<div style="text-align: center;"><figure>
<img src="layout.svg" alt="run test to generate figure" width="100%">
<figcaption>Case geometry and configuration</figcaption>
</figure></div>

The analytical solution for this case is given by[^1]

\[
    \sigma_{rr} = \frac{\sigma_y}{\sqrt{3}} \ln \left( \frac{\left( \frac{r_0}{a_0} \right)^2 + \left( \frac{a}{a_0} \right)^2 - 1}{\left( \frac{b_0}{a_0} \right)^2 + \left( \frac{a}{a_0} \right)^2 - 1} \right)
\]

where \( a_0 \) is the initial inner radius, \( b_0 \) is the initial outer radius, \( r_0 \) is the radius in the initial configuration at which the stresses are being calculated, and \( a \) is the current value of the inner radius. We consider the radial stress at the inner surface, \( r_0 = a_0 \).

The comparison between OFFBEAT and the analytical solution is shown below.

<div style="text-align: center;"><figure>
<img src="Stresses_misesPlasticity_multiMaterial.png" alt="run test to generate figure" width="75%">
<figcaption>Radial stress at the inner surface versus the inner radius of the expanding cylinder</figcaption>
</figure></div>


[^1]: E. Brunetto et al., Extension of the OFFBEAT fuel performance code to finite strains and validation against LOCA experiments, *Nuclear Engineering and Design*, Volume 406, May 2023, 112232, https://doi.org/10.1016/j.nucengdes.2023.112232