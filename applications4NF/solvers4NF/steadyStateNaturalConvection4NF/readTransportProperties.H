    Info<< "Reading transportProperties\n" << endl;

    singlePhaseTransportModel laminarTransport(U, phi);

    IOdictionary transportProperties( IOobject("transportProperties", runTime.constant(), mesh, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE));

    dimensionedScalar TRef(transportProperties.lookup("TRef"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

    // Turbulent Schmit number
    dimensionedScalar Sct(laminarTransport.lookup("Sct"));

