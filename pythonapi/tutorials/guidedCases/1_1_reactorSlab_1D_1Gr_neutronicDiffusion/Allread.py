import foamForNuclear as ffn


solver = ffn.NeutronicsSolver(region="neutroRegion", solver="SNNeutronics")

print(solver.fvSchemes.ddtSchemes)

solver.fvSchemes.import_from_openfoam()

print(solver.fvSchemes.ddtSchemes)


model2 = ffn.Model()

model2.settings.import_from_openfoam()

print(model2.settings)

model2.export_to_openfoam()
