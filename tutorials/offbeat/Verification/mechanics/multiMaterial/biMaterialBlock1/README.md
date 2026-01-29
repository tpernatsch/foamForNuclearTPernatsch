# Bi-material block 1

A 2D 10x4cm sample comprises two 5x4cm material regions subjected to a pressure of 1 MPa on the left surface. The right surface is fixed. The domain is periodic in the y-direction.


<div style="text-align: center;"><figure>
<img src="layout.svg" alt="run test to generate figure" width="100%">
<figcaption>Bi-material block case layout</figcaption>
</figure></div>

For this case the normal stress \( \sigma_{xx} \) is constant across the interface while the orthogonal stress  \( \sigma_{yy} \) is discontinuous.

<div style="text-align: center;"><figure>
<img src="sigma.png" alt="run test to generate figure" width="75%">
<figcaption>Stress distribution in the x-direction</figcaption>
</figure></div>

Conversely, the normal strain  \( \epsilon_{xx} \) is discontinuous.

<div style="text-align: center;"><figure>
<img src="epsilon.png" alt="run test to generate figure" width="75%">
<figcaption>Strain distribution in the x-direction</figcaption>
</figure></div>
