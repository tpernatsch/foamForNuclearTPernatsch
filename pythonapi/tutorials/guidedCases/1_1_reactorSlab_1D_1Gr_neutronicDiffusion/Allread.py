import foamForNuclear as ffn


solver = ffn.solvers.NeutronicsSolver(region="neutroRegion", solver="SNNeutronics")

print(solver.fvSchemes.ddtSchemes)

solver.fvSchemes.import_from_openfoam()

print(solver.fvSchemes.ddtSchemes)


model2 = ffn.case.Case()

model2.settings.import_from_openfoam()

print(model2.settings)

model2.export_to_openfoam()
