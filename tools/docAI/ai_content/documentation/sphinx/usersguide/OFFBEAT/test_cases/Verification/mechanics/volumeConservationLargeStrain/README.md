# Conservation of volume for an expanding cylinder

The [expanding cylinder case](../thickCylinderExpansion/README.md) is used to verify the conservation of volume for the large strain solver of OFFBEAT[^1]. The first two seconds of the transient are considered. The volume is plotted versus time for both the small strain and large strain solvers. The large strain solver conserves volume, while the small strain solver does not.

<div style="text-align: center;"><figure>
<img src="volumes.png" alt="run test to generate figure" width="100%">
<figcaption>Volume of the expanding cylinder versus time</figcaption>
</figure></div>

[^1]: E. Brunetto et al., Extension of the OFFBEAT fuel performance code to finite strains and validation against LOCA experiments, *Nuclear Engineering and Design*, Volume 406, May 2023, 112232, https://doi.org/10.1016/j.nucengdes.2023.112232