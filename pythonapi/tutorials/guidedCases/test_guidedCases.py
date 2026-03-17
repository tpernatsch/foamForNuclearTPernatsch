import importlib
import os


def test_1_1_reactorSlab_1D_1Gr_neutronicDiffusion(tmp_path):
    os.chdir(tmp_path)

    module = importlib.import_module("1_1_reactorSlab_1D_1Gr_neutronicDiffusion.Allrun", package=None)

    keff = module.keff

    keffRef = 1
    eps = 1e-5

    assert abs(keff - keffRef) < eps


def test_1_2_reactorSlabReflected_1D_1Gr_neutronicDiffusion(tmp_path):
    os.chdir(tmp_path)

    module = importlib.import_module("1_2_reactorSlabReflected_1D_1Gr_neutronicDiffusion.Allrun", package=None)

    keff = module.keff

    keffRef = 1
    eps = 1e-5

    assert abs(keff - keffRef) < eps


def test_1_3_reactorSphere_1D_1Gr_neutronicDiffusion(tmp_path):
    os.chdir(tmp_path)

    Allrun = importlib.import_module("1_3_reactorSphere_1D_1Gr_neutronicDiffusion.Allrun", package=None)

    # Export to OpenFOAM
    Allrun.model.export_to_openfoam()

    # Run
    Allrun.ffn.run(Allrun.model, is_preprocessing=True)

    # Post-processing
    keff = Allrun.model.keff()

    keffRef = 1
    eps = 1e-5

    assert abs(keff - keffRef) < eps


def test_1_4_reactorSlab_1D_2Gr_neutronicDiffusion(tmp_path):
    os.chdir(tmp_path)

    Allrun = importlib.import_module("1_4_reactorSlab_1D_2Gr_neutronicDiffusion.Allrun", package=None)

    # Export to OpenFOAM
    Allrun.model.export_to_openfoam()

    # Run
    Allrun.ffn.run(Allrun.model, is_preprocessing=True)

    # Post-processing
    keff = Allrun.model.keff()

    keffRef = 1
    eps = 1e-5

    assert abs(keff - keffRef) < eps


def test_1_5_reactorSlab_1D_1Gr_neutronicDiffusion_constSource(tmp_path):
    os.chdir(tmp_path)

    Allrun = importlib.import_module("1_5_reactorSlab_1D_1Gr_neutronicDiffusion_constSource.Allrun", package=None)

    # Export to OpenFOAM
    Allrun.model.export_to_openfoam()

    # Run
    Allrun.ffn.run(Allrun.model, is_preprocessing=True)

    assert True


def test_1_6_reactorSlab_1D_1Gr_neutronicDiffusion_planeSource(tmp_path):
    os.chdir(tmp_path)

    Allrun = importlib.import_module("1_6_reactorSlab_1D_1Gr_neutronicDiffusion_planeSource.Allrun", package=None)

    # Export to OpenFOAM
    Allrun.model.export_to_openfoam()

    # Run
    Allrun.ffn.run(Allrun.model, is_preprocessing=True)

    assert True


def test_2_1_reactorCylinder_2D_1Gr_neutronicDiffusion(tmp_path):
    os.chdir(tmp_path)

    Allrun = importlib.import_module("2_1_reactorCylinder_2D_1Gr_neutronicDiffusion.Allrun", package=None)

    # Export to OpenFOAM
    Allrun.model.export_to_openfoam()

    # Run
    Allrun.ffn.run(Allrun.model, is_preprocessing=True)

    # Post-processing
    keff = Allrun.model.keff()

    keffRef = 1
    eps = 1e-5

    assert abs(keff - keffRef) < eps


def test_3_1_reactorCylinder_2D_1Gr_neutronicPK_thermomech(tmp_path):
    os.chdir(tmp_path)

    Allrun = importlib.import_module("3_1_reactorCylinder_2D_1Gr_neutronicPK_thermomech.Allrun", package=None)

    # Run steady-state
    #-----------------

    # Export to OpenFOAM
    Allrun.model.export_to_openfoam()
    Allrun.ffn.run(Allrun.model, is_preprocessing=True)

    # Post-processing
    keff = Allrun.model.keff()

    keffRef = 1
    eps = 1e-5

    assert abs(keff - keffRef) < eps


    # Run transient
    #--------------

    model = Allrun.run_transient(Allrun.model)
    relErrPower = Allrun.getMaxPowerRelErr(model)

    assert abs(relErrPower) < 5
