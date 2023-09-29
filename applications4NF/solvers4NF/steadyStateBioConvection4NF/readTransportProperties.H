    Info<< "Reading transportProperties\n" << endl;

    singlePhaseTransportModel laminarTransport(U, phi);

    IOdictionary transportProperties( IOobject("transportProperties", runTime.constant(), mesh, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE));

    dimensionedScalar TRef(transportProperties.lookup("TRef"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

    // Turbulent Schmit number
    dimensionedScalar Sct(laminarTransport.lookup("Sct"));

    // Reference concentration of nanoparticle
    dimensionedScalar CRef(transportProperties.lookup("CRef"));

    // Concentration expansion coefficient
    dimensionedScalar betaC(transportProperties.lookup("betaC"));

    // Nanoparticle density
    dimensionedScalar rhop(transportProperties.lookup("rhop"));

    // Reference concentration of microorganism
    dimensionedScalar NRef(transportProperties.lookup("NRef"));

    // Concentration expansion coefficient of microorganism
    dimensionedScalar betaM(transportProperties.lookup("betaM"));

    // Microorganism density
    dimensionedScalar rhom(transportProperties.lookup("rhom"));

    // Brownian diffusion of microorganism
    dimensionedScalar Dbm(transportProperties.lookup("Dbm"));

    // Constant parameters
    dimensionedScalar W(transportProperties.lookup("W"));
    dimensionedScalar b(transportProperties.lookup("b"));

    // Constant parameters of activation energy
    dimensionedScalar kr(transportProperties.lookup("kr"));
    dimensionedScalar N(transportProperties.lookup("N"));
    dimensionedScalar Ea(transportProperties.lookup("Ea"));
