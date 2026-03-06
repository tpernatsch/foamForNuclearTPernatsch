# <b>Behavioral Models</b>

Behavioral models describe **material evolution driven by operating conditions**
such as temperature, stress, irradiation, or burnup.

Some behavioral models contribute an additional strain component (e.g. swelling
or densification), while others act through state variables or failure criteria
and do not directly modify the strain field.

In OFFBEAT, behavioral models are defined **per material** inside the
`materials` dictionary of the `solverDict` (located in the `constant` folder).

Typical examples include fuel swelling, densification, relocation, phase
transitions, and failure criteria.

???+ note
    Not all materials support all behavioral models.
    The set of available behavioral models depends on the selected material
    type, and some materials may not use any behavioral model at all.

## <b>Contents</b>

- [Densification](densification/index.md)
- [Relocation](relocation/index.md)
- [Swelling](swelling/index.md)
- [Failure Models](failure/index.md)
- [Phase Transition](phase_transition/index.md)
