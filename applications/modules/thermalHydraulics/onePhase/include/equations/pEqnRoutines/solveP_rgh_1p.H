fvScalarMatrix pEqnIncomp
(
    fvc::div(phiHbyA)
  - fvm::laplacian(rAUf, p_rgh_)
);

fvScalarMatrix pEqn(pEqnComp() + pEqnIncomp);

pEqn.setReference
(
    pRefCell_, 
    compressible ? getRefCellValue(p_rgh_, pRefCell_) : pRefValue_,
    forcePRef_
);

solverPerformance pPerf =
(
    pEqn.solve
    (
        mesh_.solver(p_rgh_.select(pimple_.finalInnerIter()))
    )
);

bool converged =
( 
    pPerf.max().nIterations() ==  0.0 
);

residual_ = max(residual_, pPerf.max().initialResidual());

if (pimple_.corrNonOrtho() == 1)
    convergedPISO = converged;

pimple_.updateFirstPISOPrevPIMPLE();