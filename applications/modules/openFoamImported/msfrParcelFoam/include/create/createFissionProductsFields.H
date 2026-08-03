// Sizes and reads the dissolved fission product species transport fields
// (Dfp_, fp_, fps_, fpp_, fpd_, gamma_, d_, qdecay_wall_). Included only
// when fpTransport_ is true, after readFissionProductsProperties.H has set
// nFPSpecies_.

{
    Info<< "\nReading fission products fields (" << nFPSpecies_ << ")" << endl;

    Dfp_.setSize(nFPSpecies_);
    fp_.setSize(nFPSpecies_);
    fps_.setSize(nFPSpecies_);
    fpp_.setSize(nFPSpecies_);
    fpd_.setSize(nFPSpecies_);
    d_.setSize(nFPSpecies_);
    qdecay_wall_.setSize(nFPSpecies_);
    gamma_.setSize(nFPSpecies_);

    forAll(fp_, i)
    {
        Dfp_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Dfp" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("Dfp" + Foam::name(i + 1), dimArea/dimTime, 0)
            )
        );

        Info<< " \n Reading field fp" << i + 1 << endl;

        fp_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fp" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        Info<< " \n Reading field fps" << i + 1 << endl;

        fps_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fp" + Foam::name(i + 1) + "s",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimless/dimVolume, 0)
            )
        );

        Info<< " \n Reading field fpp" << i + 1 << endl;

        fpp_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fp" + Foam::name(i + 1) + "p",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimless/dimVolume, 0)
            )
        );

        Info<< " \n Reading field fpd" << i + 1 << endl;

        fpd_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fp" + Foam::name(i + 1) + "d",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimless/dimArea/dimTime, 0)
            )
        );

        // Deposition velocity (output only)
        gamma_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "gamma" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("gamma" + Foam::name(i + 1), dimVelocity, 0)
            )
        );

        // Total fission products deposited
        d_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "d" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("d" + Foam::name(i + 1), dimless/dimArea, 0)
            )
        );

        // Decay power of deposited fission products
        qdecay_wall_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "qdecay_wall" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("qdecay_wall" + Foam::name(i + 1), dimPower/dimArea, 0)
            )
        );
    }
}
