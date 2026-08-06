// Reads the neutron flux / delayed-neutron precursor / decay-heat precursor
// fields, and builds the temperature/density-dependent nuclear property
// fields (D_, Sa_, Sf_, Ss_, M_), diffCoeffAlbedo_ and neutronSource_.
// Included once from the constructor body, after readNuclearProperties.H
// (and, conditionally, readFissionProductsProperties.H).

{
    Info<< "\nReading neutron flux fields (" << nEnergyGroups_ << ")" << endl;

    forAll(flux_, i)
    {
        Info<< "    Reading field flux" << i + 1 << endl;
        flux_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "flux" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        fluxTot_ += flux_[i];
    }

    Info<< "\nReading delayed neutron precursors fields (" << nPrecGroups_ << ")" << endl;

    forAll(prec_, i)
    {
        Info<< "    Reading field prec" << i + 1 << endl;
        prec_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "prec" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    Info<< "\nReading decay heat precursors fields (" << nDecGroups_ << ")" << endl;

    forAll(dec_, i)
    {
        Info<< "    Reading field dec" << i + 1 << endl;
        dec_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dec" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    Info<< "\nCreating neutron transport constants fields" << endl;

    forAll(D_, i)
    {
        word iName = Foam::name(i + 1);

        D_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "D" + iName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                D_ref_[i]
            )
        );

        Sa_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Sa" + iName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                Sa_ref_[i]
            )
        );

        Sf_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Sf" + iName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                Sf_ref_[i]
            )
        );

        forAll(D_, j)
        {
            label ijInd = nEnergyGroups_*i + j;
            word ijName = Foam::name(i + 1) + Foam::name(j + 1);

            Ss_.set
            (
                ijInd,
                new volScalarField
                (
                    IOobject
                    (
                        "Ss" + ijName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    Ss_ref_[ijInd]
                )
            );
        }

        M_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "M" + iName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                0*Nu_[i]*Sf_ref_[i]
            )
        );
    }

    forAll(neutronSource_, i)
    {
        neutronSource_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "neutronSource" + Foam::name(i + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimless/dimVolume/dimTime, 0)
            )
        );
    }
}
