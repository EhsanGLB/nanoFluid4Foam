    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


#   include "createPhi.H"
#   include "readTransportProperties.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

//************************************************ nanoFluid4Foam ************************************************//
Info<< "Reading nanoFluidProperties dictionary" << endl;
IOdictionary nanoFluidProperties(IOobject("nanoFluidProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE));

List<wordList> particlesComponents(nanoFluidProperties.lookup("particle"));
PtrList<volScalarField> alphasPtr(particlesComponents.size());

Info<< "Reading field alpha" << endl;
forAll(particlesComponents, i)
{
    alphasPtr.set(i,new volScalarField (IOobject(("alpha_"+particlesComponents[i][0]), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh));
}

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(alphasPtr, i){fields.add(alphasPtr[i]);}

Info<< "Constructing nanoFluid4Foam" << endl;
autoPtr<nanoFluidModel> nanoFluidModelPtr(nanoFluidModel::New(nanoFluidProperties, U, p, T, alphasPtr));
volScalarField rho(IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("rho", dimDensity, SMALL));
volScalarField kappa(IOobject("kappa", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("kappa", dimThermalConductivity, SMALL));
volScalarField Cp(IOobject("Cp", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL));
volScalarField mu(IOobject("mu", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("mu", dimDynamicViscosity, SMALL));
volScalarField beta(IOobject("beta", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("beta", dimThermalExpansion, SMALL));
volScalarField hm(IOobject("hm", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::NO_WRITE), mesh, dimensionedScalar("hm", dimEnthalpy, SMALL));
volScalarField Tm(IOobject("Tm", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::NO_WRITE), mesh, dimensionedScalar("Tm", dimTemperature, SMALL));
volScalarField Tmr(IOobject("Tmr", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::NO_WRITE), mesh, dimensionedScalar("Tmr", dimTemperature, SMALL));

nanoFluidModel& nanoFluidModelObj = nanoFluidModelPtr();
const PtrList<volScalarField>& DBs = nanoFluidModelObj.DBs();
const PtrList<volScalarField>& DTs = nanoFluidModelObj.DTs();

PtrList<volScalarField> moltensPtr(particlesComponents.size());
const PtrList<particleModel>& particlesProperties = nanoFluidModelObj.particlesProperties();
forAll(particlesComponents, i)
{
    moltensPtr.set(i,new volScalarField (IOobject(("molten_"+particlesComponents[i][0]), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), particlesProperties[i].molten()));
}

volScalarField kappaEff(IOobject("kappaEff", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("kappaEff", dimThermalConductivity, SMALL));
volScalarField nuEff(IOobject("nuEff", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("nuEff", dimDiffusivity, SMALL));
Info<< "Reading/calculating field rhoBB\n" << endl;
volScalarField rhoBB( IOobject( "rhoBB", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), rho*(scalar(1) - beta*(T - TRef)));

volScalarField gamma(IOobject("gamma", runTime.timeName(), mesh, IOobject::NO_READ,IOobject::AUTO_WRITE), mesh, dimensionedScalar("gamma", dimless, SMALL));

autoPtr<radiation::radiationModel> radiation(radiation::radiationModel::New(T));


