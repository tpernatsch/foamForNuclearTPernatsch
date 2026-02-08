from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile


_FVSOLUTION_SOLVER_TYPES = {
    "FPCG", "GAMG", "PBiCG", "PBiCGStab", "PCG", "PPCG", "PPCR", "smoothSolver",
    "diagonal", "MULES"
}
_FVSOLUTION_PRECONDITIONER_TYPES = {
    "DILU", "DIC", "FDIC", "GAMG", "diagonal", "distributedDIC", "none"
}
_FVSOLUTION_SMOOTHER_TYPES = {
    "DIC", "DICGaussSeidel", "FDIC", "GaussSeidel", "nonBlockingGaussSeidel",
    "symGaussSeidel"
}
_FVSOLUTION_AGGLOMERATOR_TYPES = {
    "dummy", "faceAreaPair", "algebraicPair", "assembledFaceAreaPair"
}
_FVSOLUTION_PROCESSORAGGLOMERATOR_TYPES = {"masterCoarsest"}


class fvSolutionSolver(OpenFOAMDict):
    """
    Solution solver selection.

    Parameters
    ----------
    solver : str
        Type (`FPCG`, `GAMG`, `PBiCG`, `PBiCGStab`, `PCG`, `PPCG`, `PPCR`,
        `smoothSolver`, `MULES`)
    tolerance : float
    relTol : float
    minIter : int
    maxIter : int
    preconditioner : str
        (`DILU`, `DIC`, `FDIC`, `GAMG`, `diagonal`, `distributedDIC`, `none`)
    smoother : str
        (`DIC`, `DICGaussSeidel`, `FDIC`, `GaussSeidel`, `nonBlockingGaussSeidel`, `symGaussSeidel`)
    nPreSweeps : int
    nPostSweeps : int
    nFinestSweeps : int
    cacheAgglomeration : bool
    nCellsInCoarsestLevel : int
    agglomerator : str
        (`dummy`, `faceAreaPair`, `algebraicPair`, `assembledFaceAreaPair`)
    mergeLevels : int
    processorAgglomerator : str
        (`masterCoarsest`)
    MULESCorr : bool
    nLimiterIter : int
    nAlphaCorr : int
    nAlphaSubCycles : int
    cAlpha : int
    """
    def __init__(
            self,
            solver: str,
            tolerance: float=1e-6,
            relTol: float=1e-3,
            minIter: int=None,
            maxIter: int=None,
            preconditioner: str=None,
            smoother: str=None,
            nPreSweeps: int=None,
            nPostSweeps: int=None,
            nFinestSweeps: int=None,
            cacheAgglomeration: bool=None,
            updateInterval: int=None,
            nCellsInCoarsestLevel: int=None,
            agglomerator: str=None,
            mergeLevels: int=None,
            processorAgglomerator: str=None,
            scaleCorrection: bool=None,
            directSolveCoarsest: bool=None,
            MULESCorr: bool=None,
            nLimiterIter: int=None,
            nAlphaCorr: int=None,
            nAlphaSubCycles: int=None,
            cAlpha: int=None,
            adjustSubCycles: bool=None,
            alphaMaxCo: float=None,
            solverPhase: str=None
        ):
        super().__init__()

        self.solver = solver
        self.preconditioner = preconditioner
        self.tolerance = tolerance
        self.relTol = relTol
        self.minIter = minIter
        self.maxIter = maxIter

        self.smoother = smoother
        self.nPostSweeps = nPostSweeps
        self.nPreSweeps = nPreSweeps
        self.nFinestSweeps = nFinestSweeps
        self.cacheAgglomeration = cacheAgglomeration
        self.updateInterval = updateInterval
        self.nCellsInCoarsestLevel = nCellsInCoarsestLevel
        self.agglomerator = agglomerator
        self.mergeLevels = mergeLevels
        self.processorAgglomerator = processorAgglomerator
        self.scaleCorrection = scaleCorrection
        self.directSolveCoarsest = directSolveCoarsest

        self.MULESCorr = MULESCorr
        self.nLimiterIter = nLimiterIter
        self.nAlphaCorr = nAlphaCorr
        self.nAlphaSubCycles = nAlphaSubCycles
        self.cAlpha = cAlpha
        self.adjustSubCycles = adjustSubCycles
        self.alphaMaxCo = alphaMaxCo
        self.solverPhase = solverPhase


    def __repr__(self, depth: int=0):
        if (self.solver is not None):
            super().__setitem__("solver", self.solver)
        if (self.preconditioner is not None):
            super().__setitem__("preconditioner", self.preconditioner)
        if (self.smoother is not None):
            super().__setitem__("smoother", self.smoother)
        if (self.tolerance is not None):
            super().__setitem__("tolerance", self.tolerance)
        if (self.relTol is not None):
            super().__setitem__("relTol", self.relTol)
        if (self.minIter is not None):
            super().__setitem__("minIter", self.minIter)
        if (self.maxIter is not None):
            super().__setitem__("maxIter", self.maxIter)
        if (self.nPostSweeps is not None):
            super().__setitem__("nPostSweeps", self.nPostSweeps)
        if (self.nPreSweeps is not None):
            super().__setitem__("nPreSweeps", self.nPreSweeps)
        if (self.nFinestSweeps is not None):
            super().__setitem__("nFinestSweeps", self.nFinestSweeps)
        if (self.cacheAgglomeration is not None):
            super().__setitem__("cacheAgglomeration", self.cacheAgglomeration)
        if (self.updateInterval is not None):
            super().__setitem__("updateInterval", self.updateInterval)
        if (self.nCellsInCoarsestLevel is not None):
            super().__setitem__("nCellsInCoarsestLevel", self.nCellsInCoarsestLevel)
        if (self.agglomerator is not None):
            super().__setitem__("agglomerator", self.agglomerator)
        if (self.mergeLevels is not None):
            super().__setitem__("mergeLevels", self.mergeLevels)
        if (self.processorAgglomerator is not None):
            super().__setitem__("processorAgglomerator", self.processorAgglomerator)
        if (self.scaleCorrection is not None):
            super().__setitem__("scaleCorrection", self.scaleCorrection)
        if (self.directSolveCoarsest is not None):
            super().__setitem__("directSolveCoarsest", self.directSolveCoarsest)
        if (self.MULESCorr is not None):
            super().__setitem__("MULESCorr", self.MULESCorr)
        if (self.nLimiterIter is not None):
            super().__setitem__("nLimiterIter", self.nLimiterIter)
        if (self.nAlphaCorr is not None):
            super().__setitem__("nAlphaCorr", self.nAlphaCorr)
        if (self.nAlphaSubCycles is not None):
            super().__setitem__("nAlphaSubCycles", self.nAlphaSubCycles)
        if (self.cAlpha is not None):
            super().__setitem__("cAlpha", self.cAlpha)
        if (self.adjustSubCycles is not None):
            super().__setitem__("adjustSubCycles", self.adjustSubCycles)
        if (self.alphaMaxCo is not None):
            super().__setitem__("alphaMaxCo", self.alphaMaxCo)
        if (self.solverPhase is not None):
            super().__setitem__("solverPhase", self.solverPhase)

        return(super().__repr__(depth=depth))


    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, solver):
        check_type('solver', solver, str)
        check_value('solver', solver, _FVSOLUTION_SOLVER_TYPES)
        self._solver = solver

    @property
    def preconditioner(self):
        return self._preconditioner

    @preconditioner.setter
    def preconditioner(self, preconditioner):
        check_type('preconditioner', preconditioner, (str, OpenFOAMDict), none_ok=True)
        if (isinstance(preconditioner, str)):
            check_value('preconditioner', preconditioner, _FVSOLUTION_PRECONDITIONER_TYPES)
        self._preconditioner = preconditioner

    @property
    def smoother(self):
        return self._smoother

    @smoother.setter
    def smoother(self, smoother):
        check_type('smoother', smoother, str, none_ok=True)
        if smoother is not None:
            check_value('smoother', smoother, _FVSOLUTION_SMOOTHER_TYPES)
        self._smoother = smoother

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, tolerance) -> None:
        check_type("tolerance", tolerance, (float, int), none_ok=True)
        if (tolerance is not None):
            check_positive("tolerance", tolerance)
        self._tolerance = tolerance

    @property
    def relTol(self):
        return self._relTol

    @relTol.setter
    def relTol(self, relTol) -> None:
        check_type("relTol", relTol, (float, int), none_ok=True)
        if (relTol is not None):
            check_positive("relTol", relTol)
        self._relTol = relTol

    @property
    def minIter(self):
        return self._minIter

    @minIter.setter
    def minIter(self, minIter) -> None:
        check_type("minIter", minIter, int, none_ok=True)
        if (minIter is not None):
            check_positive("minIter", minIter)
        self._minIter = minIter

    @property
    def maxIter(self):
        return self._maxIter

    @maxIter.setter
    def maxIter(self, maxIter) -> None:
        check_type("maxIter", maxIter, int, none_ok=True)
        if (maxIter is not None):
            check_positive("maxIter", maxIter)
        self._maxIter = maxIter

    @property
    def nPreSweeps(self):
        return self._nPreSweeps

    @nPreSweeps.setter
    def nPreSweeps(self, nPreSweeps) -> None:
        check_type("nPreSweeps", nPreSweeps, int, none_ok=True)
        if (nPreSweeps is not None):
            check_positive("nPreSweeps", nPreSweeps)
        self._nPreSweeps = nPreSweeps

    @property
    def nFinestSweeps(self):
        return self._nFinestSweeps

    @nFinestSweeps.setter
    def nFinestSweeps(self, nFinestSweeps) -> None:
        check_type("nFinestSweeps", nFinestSweeps, int, none_ok=True)
        if (nFinestSweeps is not None):
            check_positive("nFinestSweeps", nFinestSweeps)
        self._nFinestSweeps = nFinestSweeps

    @property
    def cacheAgglomeration(self):
        return self._cacheAgglomeration

    @cacheAgglomeration.setter
    def cacheAgglomeration(self, cacheAgglomeration) -> None:
        check_type("cacheAgglomeration", cacheAgglomeration, bool, none_ok=True)
        self._cacheAgglomeration = cacheAgglomeration

    @property
    def updateInterval(self):
        return self._updateInterval

    @updateInterval.setter
    def updateInterval(self, updateInterval) -> None:
        check_type("updateInterval", updateInterval, int, none_ok=True)
        self._updateInterval = updateInterval

    @property
    def nCellsInCoarsestLevel(self):
        return self._nCellsInCoarsestLevel

    @nCellsInCoarsestLevel.setter
    def nCellsInCoarsestLevel(self, nCellsInCoarsestLevel) -> None:
        check_type("nCellsInCoarsestLevel", nCellsInCoarsestLevel, int, none_ok=True)
        self._nCellsInCoarsestLevel = nCellsInCoarsestLevel

    @property
    def agglomerator(self):
        return self._agglomerator

    @agglomerator.setter
    def agglomerator(self, agglomerator):
        check_type('agglomerator', agglomerator, str, none_ok=True)
        if (agglomerator is not None):
            check_value('agglomerator', agglomerator, _FVSOLUTION_AGGLOMERATOR_TYPES)
        self._agglomerator = agglomerator

    @property
    def mergeLevels(self):
        return self._mergeLevels

    @mergeLevels.setter
    def mergeLevels(self, mergeLevels) -> None:
        check_type("mergeLevels", mergeLevels, int, none_ok=True)
        self._mergeLevels = mergeLevels

    @property
    def processorAgglomerator(self):
        return self._processorAgglomerator

    @processorAgglomerator.setter
    def processorAgglomerator(self, processorAgglomerator):
        check_type('processorAgglomerator', processorAgglomerator, str, none_ok=True)
        if (processorAgglomerator is not None):
            check_value('processorAgglomerator', processorAgglomerator, _FVSOLUTION_PROCESSORAGGLOMERATOR_TYPES)
        self._processorAgglomerator = processorAgglomerator

    @property
    def scaleCorrection(self):
        return self._scaleCorrection

    @scaleCorrection.setter
    def scaleCorrection(self, scaleCorrection) -> None:
        check_type("scaleCorrection", scaleCorrection, bool, none_ok=True)
        self._scaleCorrection = scaleCorrection

    @property
    def directSolveCoarsest(self):
        return self._directSolveCoarsest

    @directSolveCoarsest.setter
    def directSolveCoarsest(self, directSolveCoarsest) -> None:
        check_type("directSolveCoarsest", directSolveCoarsest, bool, none_ok=True)
        self._directSolveCoarsest = directSolveCoarsest

    @property
    def MULESCorr(self):
        return self._MULESCorr

    @MULESCorr.setter
    def MULESCorr(self, MULESCorr):
        check_type('MULESCorr', MULESCorr, bool, none_ok=True)
        self._MULESCorr = MULESCorr

    @property
    def nLimiterIter(self):
        return self._nLimiterIter

    @nLimiterIter.setter
    def nLimiterIter(self, nLimiterIter):
        check_type('nLimiterIter', nLimiterIter, int, none_ok=True)
        self._nLimiterIter = nLimiterIter

    @property
    def nAlphaCorr(self):
        return self._nAlphaCorr

    @nAlphaCorr.setter
    def nAlphaCorr(self, nAlphaCorr):
        check_type('nAlphaCorr', nAlphaCorr, int, none_ok=True)
        self._nAlphaCorr = nAlphaCorr

    @property
    def nAlphaSubCycles(self):
        return self._nAlphaSubCycles

    @nAlphaSubCycles.setter
    def nAlphaSubCycles(self, nAlphaSubCycles):
        check_type('nAlphaSubCycles', nAlphaSubCycles, int, none_ok=True)
        self._nAlphaSubCycles = nAlphaSubCycles

    @property
    def cAlpha(self):
        return self._cAlpha

    @cAlpha.setter
    def cAlpha(self, cAlpha):
        check_type('cAlpha', cAlpha, int, none_ok=True)
        self._cAlpha = cAlpha

    @property
    def adjustSubCycles(self):
        return self._adjustSubCycles

    @adjustSubCycles.setter
    def adjustSubCycles(self, adjustSubCycles):
        check_type('adjustSubCycles', adjustSubCycles, bool, none_ok=True)
        self._adjustSubCycles = adjustSubCycles

    @property
    def alphaMaxCo(self):
        return self._alphaMaxCo

    @alphaMaxCo.setter
    def alphaMaxCo(self, alphaMaxCo):
        check_type('alphaMaxCo', alphaMaxCo, (float, int), none_ok=True)
        self._alphaMaxCo = alphaMaxCo

    @property
    def solverPhase(self):
        return self._solverPhase

    @solverPhase.setter
    def solverPhase(self, solverPhase):
        check_type('solverPhase', solverPhase, str, none_ok=True)
        self._solverPhase = solverPhase


class fvSolution(OpenFOAMFile):
    """
    Solution selection.

    Parameters
    ----------
    region : str
        Name of the region

    Attributes
    ----------
    solvers : OpenFOAMDict[fvSolutionSolver]
        Dictionary containing the solvers settings.
    relaxationFactors : OpenFOAMDict
        Relation factors dict.
    extraDict : dict
        Dictionary containing additional options.
    """

    def __init__(self, region: str=""):
        super().__init__(name="fvSolution", folder="system", region=region)

        self.solvers: OpenFOAMDict[fvSolutionSolver] = OpenFOAMDict(name="solvers")

        self.relaxationFactors = OpenFOAMDict({
            "equations": OpenFOAMDict(),
            "fields": OpenFOAMDict(),
        }, name="relaxationFactors")

        self.extraDict = {}

    @property
    def relaxationFactors(self):
        return self._relaxationFactors

    @relaxationFactors.setter
    def relaxationFactors(self, relaxationFactors):
        check_type('relaxationFactors', relaxationFactors, (OpenFOAMDict, dict))

        rf = OpenFOAMDict(relaxationFactors)

        # Ensure the two subdicts are OpenFOAMDict even if user passed plain dicts
        for k in ("equations", "fields"):
            if k in rf and isinstance(rf[k], dict) and not isinstance(rf[k], OpenFOAMDict):
                rf[k] = OpenFOAMDict(rf[k])
            elif k not in rf:
                rf[k] = OpenFOAMDict()

        rf.name = "relaxationFactors"  # if you rely on name
        self._relaxationFactors = rf

    def append(self, fields: str, solver: fvSolutionSolver):
        """
        Add a new fvSolutionSolver for the specified fields.
        """
        self.add_fvSolutionSolver(fields=fields, solver=solver)


    def add_fvSolutionSolver(self, fields: str, solver: fvSolutionSolver):
        """
        Add a new fvSolutionSolver for the specified fields.
        """
        check_type("fields", fields, str)
        check_type("solver", solver, fvSolutionSolver)
        self.solvers[fields] = solver


    def add_relaxation_on_equation(self, equationName: str, value: float) -> None:
        check_type("equationName", equationName, str)
        check_type("value", value, (float, int))
        self.relaxationFactors['equations'][equationName] = value


    def add_relaxation_on_field(self, fieldName: str, value: float) -> None:
        check_type("fieldName", fieldName, str)
        check_type("value", value, (float, int))
        self.relaxationFactors['fields'][fieldName] = value


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = f"{self.solvers!r}\n"

        for key, item in self.extraDict.items():
            text += f"{key}"
            text += f"{item!r}\n"

        text += f"{self.relaxationFactors!r}\n"
        return(text)


class fvSolutionFMI(OpenFOAMFile):
    def __init__(
            self
        ):
        super().__init__("fvSolutionFMI", "system")


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter("#include", '"simulationParameters"', isAddExtraLine=True)

        text += f"libs\n(\n{tab}externalComm\n);\n\n"

        text += addParameter("host", '$host', isAddExtraLine=True)
        text += addParameter("port", '$port', isAddExtraLine=True)

        return(text)
