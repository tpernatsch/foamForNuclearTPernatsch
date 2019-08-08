Info<< "\nReading neutronics properties\n" << endl;

IOdictionary neutronicsProperties
(
    IOobject
    (
        "neutronicsProperties",
        runTime.constant(),
        neutroMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const bool eigenvalueNeutronics =
    neutronicsProperties.lookupOrDefault("eigenvalueNeutronics", false);

const bool SP3Neutronics =
    neutronicsProperties.lookupOrDefault("SP3Neutronics", false);

const bool fastNeutrons =
    neutronicsProperties.lookupOrDefault("fastNeutrons", false);

const bool adjustDiscFactors =
    neutronicsProperties.lookupOrDefault("adjustDiscFactors", false);

const bool useGivenDiscFactors =
    neutronicsProperties.lookupOrDefault("useGivenDiscFactors", false);

List<int> groupsWoDF(neutronicsProperties.lookupOrDefault<List<int> >("groupsWoDF", List<int>()));

List<int> doNotParametrize(neutronicsProperties.lookupOrDefault<List<int> >("doNotParametrize", List<int>()));
