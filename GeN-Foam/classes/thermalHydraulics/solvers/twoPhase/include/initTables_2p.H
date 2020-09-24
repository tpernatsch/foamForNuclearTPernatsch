bool writeTransferCoeffs 
(
    mesh.time().controlDict().lookupOrDefault<bool>
    (
        "writeTransferCoeffs", 
        true
    )
);

//- Drag
Kds_.insert
(
    fluid1_.name()+"."+fluid2_.name(),
    autoPtr<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                "Kd."+fluid1_.name()+"."+fluid2_.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedTensor
            (
                "", 
                dimDensity/dimTime, 
                tensor(0,0,0,0,0,0,0,0,0)
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);
Kds_.insert
(
    fluid1_.name()+".structure",
    autoPtr<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                "Kd."+fluid1_.name()+".structure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedTensor
            (
                "", 
                dimDensity/dimTime, 
                tensor(0,0,0,0,0,0,0,0,0)
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);
Kds_.insert
(
    fluid2_.name()+".structure",
    autoPtr<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                "Kd."+fluid2_.name()+".structure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedTensor
            (
                "", 
                dimDensity/dimTime, 
                tensor(0,0,0,0,0,0,0,0,0)
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);

//- Heat transfer coefficient
htcs_.insert
(
    fluid1_.name()+"."+fluid2_.name(),
    autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "htc."+fluid1_.name()+"."+fluid2_.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);

htcs_.insert
(
    fluid2_.name()+"."+fluid1_.name(),
    autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "htc."+fluid2_.name()+"."+fluid1_.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);

htcs_.insert
(
    fluid1_.name()+".structure",
    autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "htc."+fluid1_.name()+".structure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);

htcs_.insert
(
    fluid2_.name()+".structure",
    autoPtr<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "htc."+fluid2_.name()+".structure",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                (
                    (writeTransferCoeffs) ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                writeTransferCoeffs
            ),
            mesh_,
            dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);
