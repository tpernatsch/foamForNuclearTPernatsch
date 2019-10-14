{
    volScalarField& alpha1(fluid1_);
    alpha1.correctBoundaryConditions();
    volScalarField& alpha2(fluid2_);
    alpha2.correctBoundaryConditions();
    volScalarField& alphas(structure_);
    alphas.correctBoundaryConditions();
    const surfaceScalarField& phi1(fluid1_.phi());
    const surfaceScalarField& phi2(fluid2_.phi());
    word phirScheme("div(phir,alpha)");

    // Fluid1
    surfaceScalarField alphaPhic1 = 
        //fvc::flux(fluid1_.phi(), alpha1, "div(phi,alpha)"); 
        (
                fvc::flux(phi_, alpha1, "div(phi,alpha)")
            +   fvc::flux
                (
                    -fvc::flux(phi2-phi1, alpha2, phirScheme),
                    alpha1,
                    phirScheme
                )
            +   fvc::flux
                (
                    -fvc::flux(-phi1, alphas, phirScheme),
                    alpha1,
                    phirScheme
                )
        );
    fluid1_.correctInflowOutflow(alphaPhic1);
    volScalarField::Internal Sp1
    (
        IOobject
        (
            "Sp1",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Sp1", dimless/dimTime, 0)
    );
    volScalarField::Internal Su1
    (
        IOobject
        (
            "Su1",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Su1", dimless/dimTime, 0)
    );
    const volScalarField& dgdt1(fluid1_.dgdt());
    forAll(dgdt1, celli)
    {
        if (dgdt1[celli] < 0.0)
        {
            Sp1[celli] += dgdt1[celli]/max(alpha1[celli], 1e-4);
        }
        else
        {
            Su1[celli] += dgdt1[celli];
        }
    }
    myMULES::limit
    (
        geometricOneField(),
        alpha1,
        phi_,
        alphaPhic1,
        zeroField(),
        zeroField(),
        min(movingAlpha_.primitiveField(), 1.0)(),
        zeroField(),
        false
    );

    // Fluid 2
    surfaceScalarField alphaPhic2 =
        //fvc::flux(fluid2_.phi(), alpha2, "div(phi,alpha)");
        (
                fvc::flux(phi_, alpha2, "div(phi,alpha)")
            +   fvc::flux
                (
                    -fvc::flux(phi1-phi2, alpha1, phirScheme),
                    alpha2,
                    phirScheme
                )
            +   fvc::flux
                (
                    -fvc::flux(-phi2, alphas, phirScheme),
                    alpha2,
                    phirScheme
                )
        );
    fluid2_.correctInflowOutflow(alphaPhic2);
    volScalarField::Internal Sp2
    (
        IOobject
        (
            "Sp2",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Sp2", dimless/dimTime, 0)
    );
    volScalarField::Internal Su2
    (
        IOobject
        (
            "Su2",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Su2", dimless/dimTime, 0)
    );
    const volScalarField& dgdt2(fluid2_.dgdt());
    forAll(dgdt2, celli)
    {
        if (dgdt2[celli] < 0.0)
        {
            Sp2[celli] += dgdt2[celli]/max(alpha2[celli], 1e-4);
        }
        else
        {
            Su2[celli] += dgdt2[celli];
        }
    }
    myMULES::limit
    (
        geometricOneField(),
        alpha2,
        phi_,
        alphaPhic2,
        zeroField(),
        zeroField(),
        min(movingAlpha_.primitiveField(), 1.0)(),
        zeroField(),
        false
    );

    // Limit to preserve phi magnitude so that the corrected
    // alphaPhis still sum up to phi_
    {
        surfaceScalarField alpha1fu
        (
            upwind<scalar>(mesh_, phi_).interpolate(alpha1)
        );
        surfaceScalarField alpha2fu
        (
            upwind<scalar>(mesh_, phi_).interpolate(alpha2)
        );
        surfaceScalarField deltaPhi
        (
            (
                    alphaPhic1
                +   alphaPhic2
                -   phi_
            )/(alpha1fu+alpha2fu)
        );
        alphaPhic1 -= deltaPhi*alpha1fu;
        alphaPhic2 -= deltaPhi*alpha2fu;
    }

    // Solve fluid1
    fluid1_.correctInflowOutflow(alphaPhic1);
    myMULES::explicitSolve
    (
        geometricOneField(),
        alpha1,
        alphaPhic1,
        Sp1,
        Su1
    );
    alpha1 = max(alpha1, 0.0);
    fluid1_.alphaPhi() = alphaPhic1;
    fluid1_.alphaRhoPhi() = 
        fvc::interpolate(fluid1_.thermo().rho())*fluid1_.alphaPhi();

    // Solve fluid 2
    fluid2_.correctInflowOutflow(alphaPhic2);
    myMULES::explicitSolve
    (
        geometricOneField(),
        alpha2,
        alphaPhic2,
        Sp2,
        Su2
    );
    alpha2 = max(alpha2, 0.0);
    fluid2_.alphaPhi() = alphaPhic2;
    fluid2_.alphaRhoPhi() = 
        fvc::interpolate(fluid2_.thermo().rho())*fluid2_.alphaPhi();

    // Limit alphas so that they still add up to movingAlpha_
    volScalarField corr(movingAlpha_/(alpha1+alpha2));
    alpha1 *= corr;
    alpha2 *= corr;

    // Infos
    volScalarField alphaSum(alpha1 + alpha2 + alphas);
    Info<< alpha1.name() << " (avg min max) ="
        << " " << alpha1.weightedAverage(mesh_.V()).value()
        << " " << min(alpha1).value()
        << " " << max(alpha1).value()
        << endl;

    Info<< alpha2.name() << " (avg min max) ="
        << " " << alpha2.weightedAverage(mesh_.V()).value()
        << " " << min(alpha2).value()
        << " " << max(alpha2).value()
        << endl;

    Info<< "sum(alpha)" << " (avg min max) ="
        << " " << alphaSum.weightedAverage(mesh_.V()).value()
        << " " << min(alphaSum).value()
        << " " << max(alphaSum).value()
        << endl;
}
