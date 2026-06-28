# **Multi-Material Interface**

While further details will be covered in upcoming publications and in the Theory Manual, the multi-material interface was originally developed to enable the treatment of multi-material bodies with perfectly bonded interfaces (i.e., no external patch separating different `cellZones`) in OFFBEAT.

This treatment is inspired by the original work of Tukovic et al. [@TukovicMultiMaterial] (as well as by the implicit contact boundary condition developed by Scolaro et al. [@ScolaroImplicitContact]). It extends their approach within the mesh.

It has been observed that the `multiMaterialInterface` correction also enhances accuracy in single-material domains, particularly in cases with significant temperature gradients, variable material properties, or additional strain components such as creep. Indeed, when force conservation is applied fully between adjacent cells with such non-uniformities, it results in a different stress discretization (for both the implicit and explicit components), as well as a different internal face displacement interpolation (and thus a different gradient scheme). **This correction is not activated by default, but it is recommended to activate it in most cases, if convergence permits.**

In brief, without delving into the mathematical details, the `multiMaterialCorrection`, when used in a `mechanicsSubSolver`, employs custom gradient and Laplacian schemes that enforce full conservation at cell interfaces. **Consequently, the schemes specified in `fvSchemes` for the gradient and Laplacian for the displacement fields are disregarded in favor of this corrected approach.**

The `multiMaterialInterface` correction is activated automatically when the `multiMaterialCorrection` dictionary is defined in the `mechanicsOptions` dictionary, which in turn is located within the `solverDict`. The type of the `multiMaterialInterface` class must then be selected with the `<code><b>type</b></code>` keyword.

## <b>Classes</b>

???+ warning

    Although the reasons are not entirely understood, it has been found that applying the full multi-material correction for every internal face may negatively impact convergence. For this reason, the `multiMaterialInterface` classes use a system of weights to allow for a partial use of the correction.

Currently, OFFBEAT supports the following mechanics solvers:

- [uniform](../../classes/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniformMultiMaterialInterface.md) - It applies the multi-material correction uniformly in the mesh (a good choice for standard fuel performance simulations where there are no multi-material bodies).
- [uniform2D](../../classes/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniform2DMultiMaterialInterface.md) - It applies the multi-material correction uniformly in the mesh, but only in the x-y plane (a good choice for standard fuel performance simulations that are quite coarse along the axial direction, where there are no multi-material bodies).
- [cellZone](../../classes/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/cellZoneMultiMaterialInterface.md) - The user can select `defaultWeights` for the correction on internal faces (within a cellZone) and an `interfaceWeight` for the interface faces, defined as the faces between two different cellZones.
- [faceSet](../../classes/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/faceSetMultiMaterialInterface.md) - The user can select `defaultWeights` for the correction on internal faces (within a cellZone) and an `interfaceWeight` for the interface faces, defined using `faceSet`.
- [topoSet](../../classes/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/topoSetSourceMultiMaterialInterface.md) - The user can select `defaultWeights` for the correction to be applied to internal faces (within a cellZone) and an `interfaceWeight` for the interface faces, defined using `topoSet`.