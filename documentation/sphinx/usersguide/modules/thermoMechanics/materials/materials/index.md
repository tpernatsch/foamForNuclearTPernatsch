# <b>Materials</b>

In OFFBEAT, a **material** defines the physical characteristics of the solid
assigned to a given mesh region.

Materials are associated with mesh regions through `cellZone`s and act as
containers for:
- material properties (thermal and mechanical),
- behavioral models (e.g. swelling, densification),
- constitutive behavior (e.g. elasticity, plasticity),
- optional damage and failure models.

Materials are selected using the <code><b>type</b></code> keyword inside each
material subdictionary of the <code><b>materials</b></code> section in the
<code><b>solverDict</b></code>.

Each material comes with a **default set of thermo-mechanical properties and
associated models**, which reflect typical correlations and assumptions for that
material (e.g. UO₂, Zircaloy, steels).

These defaults can be selectively overridden by the user through explicit model
choices in the `materials` dictionary.

---

## Material selection

```cpp
materials
{
    fuel
    {
        type UO2;
        // ...
        // material properties, behavior, constitutive models, etc.
    }
}
```

The material type defines the **base material class** (e.g. `UO2`, `Zircaloy`,
`Steel1515Ti`), independently of the `cellZone` name.

---

## Scope of this section

This section documents the **available material types** in OFFBEAT and their
intended usage.

Specific details on the various models for material properties, behavioral phenomena, constitutive behavior, damage and failure, are documented in the corresponding subsections of the user guide.

---

## <b>Classes</b>

Currently, OFFBEAT supports the following material types:

* [constant](../../../classes/materials/materialModel/constantMaterial.md)
  generic material with user-defined constant properties.

* [UO2](../../../classes/materials/materialModel/UO2.md)
  uranium dioxide fuel material.

* [UPuO2](../../../classes/materials/materialModel/UPuO2.md)
  mixed oxide (MOX) fuel material.

* [Zircaloy](../../../classes/materials/materialModel/Zircaloy.md)
  Zircaloy cladding material.

* [Steel1515Ti](../../../classes/materials/materialModel/Steel1515Ti.md)
  15-15Ti austenitic steel.

* [HastelloyN](../../../classes/materials/materialModel/HastelloyN.md)
  Hastelloy N alloy.

* [Inconel600](../../../classes/materials/materialModel/Inconel600.md)
  Inconel 600 alloy.

* [Molybdenum](../../../classes/materials/materialModel/Molybdenum.md)
  molybdenum structural material.

* [SiC](../../../classes/materials/materialModel/SiC.md)
  silicon carbide material.

* [PyC](../../../classes/materials/materialModel/PyC.md)
  pyrolytic carbon material.

* [Buffer](../../../classes/materials/materialModel/Buffer.md)
  buffer material for coated particle fuels.