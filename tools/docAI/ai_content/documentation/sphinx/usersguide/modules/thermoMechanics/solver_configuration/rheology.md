# <b>Rheology model</b>

The <code><b>rheology</b></code> class is designed to handle the material's rheological properties, or constitutive mechanical model, in OFFBEAT.

With the use of simple switches, it also allows the user to activate or deactivate thermal strains, the plane stress approximation (useful, for example, for r-theta simulations), or the modified plane strain approximation (useful for 1.5D simulations).

## <b>Classes</b>

Currently, OFFBEAT supports only the following option, which is **the default choice for OFFBEAT**:

- [byMaterial](../../classes/rheology/rheologyByMaterial.md) - The constitutive laws for each material (e.g. elasticity, plasticity, etc.) are selected at the material level. This is done by choosing the appropriate value for the keyword `rheologyModel` in each `materials` subdictionary.

<!-- The main constitutive laws (to be selected one for each material) currently available in OFFBEAT are:

- [elasticity](@ref Foam.elasticity), linear elasticity following Hooke's law is assumed. *Note that even if the constitutive law is linear-elastic, the total strain field might include additional components, such as thermal strain or swelling.*
- [misesPlasticity](@ref Foam.misesPlasticity), following Von Mises theory an additional (instantaneous) plastic strain component is calculated as a function of the stress state and of the cumulated plastic strain itself. Hooke's law is assumed to calculate the stress as a function of the elastic strain component.
- [misesPlasticCreep](@ref Foam.misesPlasticCreep), an additional time-dependent plastic strain component is calculated according to the specific creep model selected. Von Mises theory is used for the instantaneous portion of the plastic strain, while Hooke's law is assumed for calculating the stress as a function of the (remaining) elastic strain component.

**The previous models can be used both for linear and non-linear mechanical solvers (i.e. they are suitable both for the small strain approximation and for the finite-strain framework).**
The following ones instead should be used only in combination with large-strain mechanical solvers:

- [hyperElasticity](@ref Foam.hyperElasticity), the stress follows a StVenant-Kirchhoff hyperElasticity constitutive law.
- [hyperElasticMisesPlasticCreep](@ref Foam.hyperElasticMisesPlasticCreep), a time-dependent plastic strain component is calculated according to the specific creep model selected. The total strain is considered to be additive. The creep and additional strain component are removed from the total strain tensor, and the stress is calculated following a StVenant-Kirchhoff hyperElasticity constitutive law as a function of the (remaining) elastic strains.
- [neoHookeanElasticity](@ref Foam.neoHookeanElasticity), the stress follows a neoHookean elasticity constitutive law.
- [neoHookeanMisesPlasticity](@ref Foam.neoHookeanMisesPlasticity), Mises plasticity is used for the plastic component, while the stress follows a neoHookean elasticity constitutive law.


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#cfe3f5; color:#05134a'>
<b>Note on the correct laws to use for large-strain simulations</b>
br><br>When modeling a fuel rod in a large strain framework, it is typically still a good approximation to select the `misesPlasticCreep` law if one wants to include plastic phenomena in the simulation (e.g. for the cladding).
<br><br>Neo-Hookean or hyperelastic laws are strongly not suggested (they do not represent the behavior of fuel and cladding). They are present in OFFBEAT mainly for benchmark and verification purposes.
</div> -->