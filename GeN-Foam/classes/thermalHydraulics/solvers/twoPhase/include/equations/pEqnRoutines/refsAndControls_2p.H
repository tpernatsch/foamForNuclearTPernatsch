surfaceScalarField& phi1(fluid1_.phi());
surfaceScalarField& phi2(fluid2_.phi());
dimensionedScalar residualAlpha1(fluid1_.residualAlpha());
dimensionedScalar residualAlpha2(fluid2_.residualAlpha());
surfaceScalarField alphaf1(fvc::interpolate(alpha1));
surfaceScalarField alphaf2(fvc::interpolate(alpha2));
surfaceScalarField movingAlphaf(fvc::interpolate(movingAlpha_));

bool noPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::none
);
bool explicitPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::legacy
);
bool implicitPartialElimination
(
    partialEliminationMode_ == partialEliminationMode::implicit
);
