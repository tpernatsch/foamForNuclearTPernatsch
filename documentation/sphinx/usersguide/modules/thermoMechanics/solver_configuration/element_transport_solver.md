# Elements Transport Solvers {#elementTransport}

OFFBEAT features transport solvers for the solution of elements redistribution within the fuel rod domain. If one (or more) of these solvers are enabled by the user, the solution of one (or more) PDEs will be included in the OFFBEAT solution loop. 

The element transport solver <u>must</u> be selected with the <code><b>elementTransport</b></code> keyword in the main dictionary of OFFBEAT (i.e. the <code><b>solverDict</b></code>  dictionary, located in the <code><b>constant</b></code>  folder). 
<br>

The two following classes allow one to handle element transport solvers in the OFFBEAT simulation:

- [fromLatestTime](@ref Foam.elementTransport), that disables all the element transport solvers.
- [byList](@ref Foam.elementTransportByList), that allows one to specify a list of element transport solvers that the user wants to include in the fuel performance simulation. The following transport solvers are currently implemented in OFFBEAT:

	- [porosityTransport](@ref Foam.porosityTransportSolver), for the porosity migration in nuclear fuel pins.
	- [PuRedistribution](@ref Foam.PuRedistributionSolver), for the Plutonium redistribution during irradiation in fuel pins.
	- [AmRedistribution](@ref Foam.AmRedistributionSolver), for the Americium redistribution during irradiation in fuel pins.

***
## The porosity migration solver

The radial redistribution of porosity in fast-MOX fuels is responsible for an early restructuring of the pellets taking place already from the first hours of irradiation. In particular, porosity migration induces the formation and growth of a central void, surrounded by a region where the trails left by the migrating pores cause the formation of a columnar grain region. As fuel restructuring determines local changes in thermo-mechanical properties of the fuel, an accurate assessment of the extensions of these zones is crucial for fuel behavior analyses. The governing equation for porosity migration implemented in OFFBEAT is:

\[
\frac{\partial p}{\partial t} + \nabla \cdot\left(p\left(1-p\right) \vec{v_p}\right) - \nu\nabla^2p= 0
\tag{1}
\]

The fuel fractional porosity is represented by the scalar quantity \( p \), which is transported through the domain by an advective contribution (i.e., the divergence term) and a diffusive contribution (i.e., the Laplacian term). The advective term is responsible for modeling the pores' evaporation-condensation effect through the pore velocity \( \vec{v_p} \), which is provided by specific correlations such as the one by Sens for UO\(_2\) fuels (Sens, 1972):

\[
\vec{v_p} = a 
\left(
    b + 
    c~T + 
    d~T^2 + 
    e~T^3
\right)
T^{-2.5}\Delta H p_0~\exp(-\Delta H/RT)
~\vec{\nabla T}
\tag{2}
\]

where \( a,b,c,d \) and \( e \) are model constants, \( R \) is the universal gas constant, and \( p_0 \) and \( \Delta H \) are the vapor pressure and the enthalpy of vaporization. Besides the Sens model, other pore velocity correlations implemented in OFFBEAT are the Lackey (Lackey, 1972) or the Clement-Finnis model (Clement, 1978) for MOX fuels and the more advanced approach proposed by Ikusawa (Ikusawa, 2014), where MOX pore velocity is determined accounting for the vapor pressures of different fuel constituents.

With the local modification to fuel thermal properties (namely thermal conductivity and volumetric heat source) with porosity evolution, pore velocity locally tends to zero for the zones characterized by a porosity fraction that approaches 1. Despite the null velocity ensuring a physical limit to the porosity fraction, this bound is further enforced by the presence of the \( (1 - p) \) term multiplied by \( p \) in the divergence of Eq. (1). Besides limiting the maximum porosity to 100%, this term also contributes to reducing oscillations in the porosity transport solution. With a similar intent, the \( -\nu\nabla^2 p \) term is added to Eq. (1), typically considering diffusivities \( \nu \) in the order of \( 10^{-12} \) to \( 10^{-10} \) m²/s. Despite having a minor contribution to porosity redistribution, the diffusive term is included to account for surface and bulk diffusion, as stated by the authors of (Novascone et al., 2018).

To facilitate convergence and avoid numerical oscillations in the solution of the porosity migration equation, OFFBEAT employs a time-stepping criterion based on the Courant-Friedrichs-Lewy (CFL) condition, commonly used for solving advective partial differential equations, particularly in fluid dynamics. The condition, shown in Eq. (3), ensures that within a time step \( \Delta t \), the transported information (in our case, the scalar quantity \( p \) carried by the velocity \( \vec{v_p} \)) does not travel farther than the grid cell size \( \Delta x \). This prevents solution instability and non-physical behavior. Once the user-specified value for the maximum Courant number \( C_{\mathrm{max}} \) is known (typically a limit around 0.5 is used), OFFBEAT computes the maximum allowable \( \Delta t \) by inverting Eq. (3).

\[
\frac{v_p \Delta t}{\Delta x} \leq C_{max}
\tag{3}
\]

### Usage

The porosity migration solver can be enabled from `constant/solverDict`:
```cpp

elementTransport byList;

elementTransportOptions
{
    solvers
	( 
		porosityTransport 
	);

	porosityOptions
	{
		diffusiveTerm   		off; // Default: off
		boundedAdvectionTerm	off; // Default: off
	}
}
```
The maximum allowed courant number can be specified as a time-stepping criterion from `system/controlDict` as:
```cpp
maxCourantNoPorosity	0.5;

```

***
## Plutonium and Americium redistribution solvers

Actinide redistribution in MOX fuel involves the migration of heavy metal atoms, such as plutonium, americium, and neptunium, up the temperature gradient towards hotter regions of the fuel pellet. The migration of actinides is driven by two concurrent phenomena: **thermodiffusion** and **vapor transport**. Thermodiffusion, also known as the **Soret effect**, is triggered in the presence of a temperature gradient and refers to the migration of heavier particles toward colder regions. The mass flux caused by the temperature gradient can be modeled using the equation shown in Eq. (1). In Eq. (1), the flux \( \vec{J}_{\text{Diff}} \) has two components: a Fickian contribution driven by the concentration gradient, and a Soret contribution driven by the temperature gradient. In this equation, \( c \) represents the mass fraction of a generic actinide, \( D \) is the thermal diffusion coefficient, \( Q \) is the molar heat of transport, \( R \) is the universal gas constant, and \( T \) is the temperature.

\[
\vec{J}_{\text{Diff}} = -D\left(\vec{\nabla c} + c\left(1-c\right)\frac{Q}{RT^2}\vec{\nabla T}\right)
\tag{1}
\]

The other contribution to actinide migration arises from the transport of gaseous species due to migrating porosity. **Vapor transport** is significant during the initial phase of irradiation, when fuel restructuring occurs, but thermodiffusion becomes the dominant mechanism over long-term irradiation. The mass flux due to pore migration is represented by Eq. (2), where \( D \) is the diffusion coefficient, \( A \) is a constant, \( p \) is the porosity, \( l \) and \( d \) are the characteristic pore diameter and thickness, respectively, and \( v_p \) is the pore velocity.

\[
\vec{J}_{\text{vapor}} = -D\cdot A p\frac{l}{d}\vec{\nabla T} \exp\left(-\frac{D}{l v_p}\right)\cdot c
\tag{2}
\]

Further details on the mechanisms governing these two migration contributions can be found in Olander's manual (Olander, 1976), while information on the determination of the constants involved can be found in (Di Marcello et al., 2014).

The spatial distribution of the mass concentration of a specific actinide, \( c \), is calculated in **OFFBEAT** for each timestep by solving the partial differential equation shown in Eq. (3):

\[
\frac{\partial c}{\partial t} + \nabla \cdot\left(\vec{J}_{\text{Diff}}+\vec{J}_{\text{vapor}}\right) = 0
\tag{3}
\]

To ensure the conservation of actinides' mass within the fuel domain, a zero-flux boundary condition must be applied at the domain boundaries when solving Eq. (3). Considering the typical radial temperature profile in a nuclear fuel pellet, the vapor transport outflow at the fuel surface is expected to be zero, which naturally enforces this condition. This is because the temperature at the outer surface of the fuel is sufficiently low to prevent porosity migration, and even in the presence of annular pellets, where an inner fuel surface exists, the temperature gradient (and hence the pore velocity) is expected to be negligible or zero.

Regarding **Soret diffusion**, due to the presence of a significant temperature gradient at the fuel's outer surface, a **zero-current constraint** must be enforced. To achieve this, a new boundary condition has been developed in **OFFBEAT**. This boundary condition, called `zeroCurrentActinidesRedistributionFvPatchScalarField`, has been derived from the OpenFOAM class `fixedGradientFvPatchScalarField`, which allows for imposing the gradient of a scalar quantity at a specific domain boundary.

Considering a generic boundary with normal vector \( \vec{n} \), by imposing \( \vec{J}_{\text{Diff}} \cdot \vec{n} = 0 \), the surface normal concentration gradient at the domain boundary \( \vec{\nabla c}\cdot \vec{n} \) can be expressed as:

\[
\vec{\nabla c} \cdot \vec{n} = c~(c - 1)~\frac{Q}{RT^2}\vec{\nabla T} \cdot \vec{n} 
\tag{4}
\]


### Usage

The actinides redistribution solvers can be enabled from `constant/solverDict`:
```cpp

elementTransport byList;

elementTransportOptions
{
    solvers
    (
        AmRedistribution
        PuRedistribution
    );

    "AmRedistributionOptions|PuRedistributionOptions"
    {
        poreMigration off;
    }
}
```

***
## References

- Sens, P. F. (1972). *The kinetics of pore movement in UO₂ fuel rods*. Journal of Nuclear Materials, **43**(3), 293–307. [https://doi.org/10.1016/0022-3115(72)90061-X](https://doi.org/10.1016/0022-3115(72)90061-X)

- Lackey, W. J., Homan, F. J., & Olsen, A. R. (1972). *Porosity and actinide redistribution during irradiation of (U,Pu)O₂*. Oak Ridge National Laboratory (ORNL), Oak Ridge, TN. [https://doi.org/10.2172/4654410](https://doi.org/10.2172/4654410)

- Clement, C. F., & Finnis, M. W. (1978). *The movement of lenticular pores in mixed oxide (U, Pu)O₂ nuclear fuel elements*. Journal of Nuclear Materials, **75**(1), 115–124. [https://doi.org/10.1016/0022-3115(78)90035-1](https://doi.org/10.1016/0022-3115(78)90035-1)

- Ikusawa, Y., Ozawa, T., Hirooka, S., Maeda, K., Kato, M., & Maeda, S. (2014). *Development and verification of the thermal behavior analysis code for MA containing MOX fuels*. In *Proceedings of the 22nd International Conference on Nuclear Engineering (ICONE22)*, Volume 1: Plant Operations, Maintenance, Engineering, Modifications, Life Cycle and Balance of Plant; Nuclear Fuel and Materials; Plant Systems, Structures and Components; Codes, Standards, Licensing and Regulatory Issues. American Society of Mechanical Engineers. [https://doi.org/10.1115/ICONE22-30005](https://doi.org/10.1115/ICONE22-30005)

- Novascone, S., Medvedev, P., Peterson, J. W., Zhang, Y., & Hales, J. (2018). *Modeling porosity migration in LWR and fast reactor MOX fuel using the finite element method*. Journal of Nuclear Materials, **508**, 226–236. [https://doi.org/10.1016/J.JNUCMAT.2018.05.041](https://doi.org/10.1016/J.JNUCMAT.2018.05.041)

### References

- Meyer, R. O. (1974). *Analysis of plutonium segregation and central-void formation in mixed-oxide fuels*. Journal of Nuclear Materials, **50**(1), 11–24. [https://doi.org/10.1016/0022-3115(74)90055-5](https://doi.org/10.1016/0022-3115(74)90055-5)

- Olander, D. R. (1976). *Fundamental Aspects of Nuclear Reactor Fuel Elements*. Technical Information Center of Public Affairs. ISBN: 0-87079-031-5.

- Di Marcello, V., Rondinella, V., Schubert, A., Van De Laar, J., & Van Uffelen, P. (2014). *Modelling actinide redistribution in mixed oxide fuel for sodium fast reactors*. Progress in Nuclear Energy, **72**, 83–90. Elsevier.

- Di Marcello, V., Schubert, A., Van De Laar, J., & Van Uffelen, P. (2012). *Extension of the TRANSURANUS plutonium redistribution model for fast reactor performance analysis*. Nuclear Engineering and Design, **248**, 149–155. [https://doi.org/10.1016/J.NUCENGDES.2012.03.037](https://doi.org/10.1016/J.NUCENGDES.2012.03.037)

***
Return to [User Manual](@ref userManual)
