# Force conservation

We consider the traditional Hertz contact problem of a 5 cm radius cylinder pressing downward on a flat surface[^1], as illustrated below. The purpose of this case is not to verify the correct solution of the Hertz problem, but rather to demonstrate force conservation for a case in which the surfaces are not aligned.

The mechanical properties of the flat base and cylinder are \(E=50\) GPa and \(\nu=0.3\).

<div style="text-align: center;"><figure>
<img src="layout.svg" alt="run test to generate figure">
<figcaption>Case geometry and configuration</figcaption>
</figure></div>

The analytical solution for the contact width \(b\), the maximum contact pressure \(P_{max}\), and the stress distribution in the \(y\)-/vertical direction \(\sigma_{yy}(y)\), from [^1], is

\[
    b = \sqrt{\frac{8}{\pi} F''R^2 \frac{1-\nu^2}{E}}
\]

\[
    P_{max} = \frac{2F''R}{\pi b}
\]

\[
    \sigma_{yy}(y) = \frac{-P_{max}}{\sqrt{1+\frac{y^2}{b^2}}}
\]

The comparison between the OFFBEAT solution and the analytical solution is shown in the figure below.

<div style="text-align: center;"><figure>
<img src="sigma.png" alt="run test to generate figure">
<figcaption>Comparison between OFFBEAT and analytical solution for the Hertz contact problem</figcaption>
</figure></div>

[^1]: R. G. Budynas and J. K. Nisbett, "Shigley's Mechanical Engineering Design, 11th edition," McGraw-Hill, Jan 2019.