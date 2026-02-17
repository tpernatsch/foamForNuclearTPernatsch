import foamForNuclear as ffn

ffn.allclean()

from model import *

print(model)

model.export_to_openfoam()

#==============================================================================*

# ffn.copyFolder('dynamicMeshDict', 'constant')

ffn.run_preprocessing(model)

model.plot_mesh(region=thMesh, show_edges=True)
model.plot_mesh(region=thMesh, show_edges=True, fieldName='alpha.water')

ffn.run(model)

#==============================================================================*
# Post-processing

from Allpostprocess import *

#==============================================================================*
