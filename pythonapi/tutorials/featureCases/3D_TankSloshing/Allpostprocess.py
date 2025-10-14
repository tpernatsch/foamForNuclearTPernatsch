from model import *

try:
    model.plot_mesh(
        region=thMesh,
        time=model.settings.endTime,
        show_edges=True,
        fieldName='alpha.water'
    )
except:
    pass

model.plot_animation(
    region="",
    fieldName="U",
    cmap='Blues',
    threshold=[0.9, 1.1],
    thresholdFieldName="alpha.water",
    unit='m/s',
    fps=5
)