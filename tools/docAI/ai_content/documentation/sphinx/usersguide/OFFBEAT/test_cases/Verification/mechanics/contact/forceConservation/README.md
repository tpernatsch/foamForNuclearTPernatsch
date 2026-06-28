# Force conservation

We consider the traditional Hertz contact problem of a 5 cm radius cylinder pressing downward
onto a flat surface[^1], as illustrated below. The purpose of this case is not to verify the
correct Hertz solution. Instead, it is to demonstrate force conservation for a situation in
which the surfaces are not aligned.

The mechanical properties of the flat base are \(E=50\) GPa and \(\nu=0.3\). For the cylinder,
we use \(E=200\) GPa and \(\nu=0.3\).

For the OFFBEAT analysis, we deliberately choose mismatching nodalisations at the contact
point. As a result, corrections are required to ensure conservation of forces on the
contact surfaces. The top, bottom, and contact surfaces should each have a total force
equal to the top pressure multiplied by the surface area of the top surface. The force
must act in the vertical direction.

<div style="text-align: center;"><figure>
<img src="layout.svg" alt="run test to generate figure">
<figcaption>Case geometry and configuration</figcaption>
</figure></div>

[^1]: R. G. Budynas and J. K. Nisbett, "Shigley's Mechanical Engineering Design, 11th edition," McGraw-Hill, Jan 2019.