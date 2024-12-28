    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
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

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "alpha1");

//------------------- NF4Foam ------------------//
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

    //------------ phase1 ------------//
    Info<< "Reading nanoFluidProperties1 dictionary" << endl;
    IOdictionary nanoFluidProperties1
    (
        IOobject
        (
            "nanoFluidProperties1",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<wordList> particlesComponents1(nanoFluidProperties1.lookup("particle"));
    PtrList<volScalarField> alphasPtr1(particlesComponents1.size());

    Info<< "Reading field alphap" << endl;
    forAll(particlesComponents1, i)
    {
        alphasPtr1.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    ("alpha1_"+particlesComponents1[i][0]),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields1;
    forAll(alphasPtr1, i)
    {
        fields1.add(alphasPtr1[i]);
    }

    Info<< "Constructing nanoFluid4Foam" << endl;
    autoPtr<nanoFluidModel> nanoFluidModelPtr1(nanoFluidModel::New(nanoFluidProperties1, U, pd, T, alphasPtr1));

    Info<< "Constructing field rho1\n" << endl;
    volScalarField rho1
    (
        IOobject
        (
            "rho1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho1", dimDensity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field Cp1\n" << endl;
    volScalarField Cp1
    (
        IOobject
        (
            "Cp1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp1", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field kappa1\n" << endl;
    volScalarField kappa1
    (
        IOobject
        (
            "kappa1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa1", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field mu1\n" << endl;
    volScalarField mu1
    (
        IOobject
        (
            "mu1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu1", dimDynamicViscosity, SMALL),
        "zeroGradient"
    );

    nanoFluidModel& nanoFluidModelObj1 = nanoFluidModelPtr1();
    const PtrList<volScalarField>& DBs1 = nanoFluidModelObj1.DBs();
    const PtrList<volScalarField>& DTs1 = nanoFluidModelObj1.DTs();

    PtrList<volScalarField> moltensPtr1(particlesComponents1.size());
    const PtrList<particleModel>& particlesProperties1 = nanoFluidModelObj1.particlesProperties();
    forAll(particlesComponents1, i)
    {
        moltensPtr1.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    ("molten1_"+particlesComponents1[i][0]),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                particlesProperties1[i].molten()
            )
        );
    }

    //------------ phase2 ------------//
    Info<< "Reading nanoFluidProperties2 dictionary" << endl;
    IOdictionary nanoFluidProperties2
    (
        IOobject
        (
            "nanoFluidProperties2",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<wordList> particlesComponents2(nanoFluidProperties2.lookup("particle"));
    PtrList<volScalarField> alphasPtr2(particlesComponents2.size());

    Info<< "Reading field alphanp" << endl;
    forAll(particlesComponents2, i)
    {
        alphasPtr2.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    ("alpha2_"+particlesComponents2[i][0]),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields2;
    forAll(alphasPtr2, i)
    {
        fields2.add(alphasPtr2[i]);
    }

    Info<< "Constructing nanoFluid4Foam" << endl;
    autoPtr<nanoFluidModel> nanoFluidModelPtr2(nanoFluidModel::New(nanoFluidProperties2, U, pd, T, alphasPtr2));

    Info<< "Constructing field rho2\n" << endl;
    volScalarField rho2
    (
        IOobject
        (
            "rho2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho2", dimDensity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field Cp2\n" << endl;
    volScalarField Cp2
    (
        IOobject
        (
            "Cp2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp2", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field kappa2\n" << endl;
    volScalarField kappa2
    (
        IOobject
        (
            "kappa2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa2", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    Info<< "Constructing field mu2\n" << endl;
    volScalarField mu2
    (
        IOobject
        (
            "mu2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu2", dimDynamicViscosity, SMALL),
        "zeroGradient"
    );

    nanoFluidModel& nanoFluidModelObj2 = nanoFluidModelPtr2();
    const PtrList<volScalarField>& DBs2 = nanoFluidModelObj2.DBs();
    const PtrList<volScalarField>& DTs2 = nanoFluidModelObj2.DTs();

    PtrList<volScalarField> moltensPtr2(particlesComponents2.size());
    const PtrList<particleModel>& particlesProperties2 = nanoFluidModelObj2.particlesProperties();
    forAll(particlesComponents2, i)
    {
        moltensPtr2.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    ("molten2_"+particlesComponents2[i][0]),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                particlesProperties2[i].molten()
            )
        );
    }

//------------------- interFoam ------------------//

    rho1 = nanoFluidModelPtr1->rho();
    rho2 = nanoFluidModelPtr2->rho();
    Cp1 = nanoFluidModelPtr1->Cp();
    Cp2 = nanoFluidModelPtr2->Cp();

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho1)*phi
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + rho*gh
    );


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, pimple.dict(), pdRefCell, pdRefValue);
    mesh.schemesDict().setFluxRequired(pd.name());

    scalar pRefValue = 0.0;

    if (pd.needReference())
    {
        pRefValue = readScalar(pimple.dict().lookup("pRefValue"));

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pdRefCell)
        );
    }

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

//------------------- NF4Foam -> 2phase ------------------//
    /*Info<< "Reading field Cp\n" << endl;
    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cp", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );*/

    /*Info<< "Reading field kappa\n" << endl;
    volScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimThermalConductivity, SMALL),
        "zeroGradient"
    );*/

    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1*Cp1 + (scalar(1) - alpha1)*rho2*Cp2,
        alpha1.boundaryField().types()
    );
    rhoCp.oldTime();

    surfaceScalarField rhoCpPhi
    (
        IOobject
        (
            "rhoCpPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rhoCp)*phi
    );

//------------------- NF4Foam -> Turbulence ------------------//
    Info<< "Constructing field kappat\n" << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappat", dimThermalConductivity, SMALL)
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    dimensionedScalar Prt(twoPhaseProperties.lookup("Prt"));
    dimensionedScalar Sct(twoPhaseProperties.lookup("Sct"));


//------------------- NF4Foam -> Porous Media ------------------//
    porousZones pZones(mesh);

    IOdictionary porousZonesDict
    (
        IOobject
        (
            "porousZones",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("porosity", dimless, 1.0),
        "zeroGradient"
    );

    volScalarField rhoPM
    (
        IOobject
        (
            "rhoPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoPM", dimDensity, SMALL),
        "zeroGradient"
    );

    volScalarField CpPM
    (
        IOobject
        (
            "CpPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("CpPM", dimSpecificHeatCapacity, SMALL),
        "zeroGradient"
    );

    volScalarField kappaPM
    (
        IOobject
        (
            "kappaPM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaPM", dimThermalConductivity, SMALL),
        "zeroGradient"
    );

    forAll( pZones, zoneI )
    {
        const label & zoneId( pZones[zoneI].zoneId() );
        const labelList & cells(mesh.cellZones()[zoneId]);
        const scalar & zonePorosity( pZones[zoneI].porosity() );

        const dictionary& thermalPropertiesPM = pZones[zoneI].dict().subDict("thermalProperties");
        dimensionedScalar zoneRho(thermalPropertiesPM.lookup("rho"));
        dimensionedScalar zoneCp(thermalPropertiesPM.lookup("Cp"));
        dimensionedScalar zoneKappa(thermalPropertiesPM.lookup("kappa"));

        forAll( cells, cellI )
        {
            porosity[cells[cellI]] = zonePorosity;
            rhoPM[cells[cellI]] = zoneRho.value();
            CpPM[cells[cellI]] = zoneCp.value();
            kappaPM[cells[cellI]] = zoneKappa.value();
        }
    }

//------------------- NF4Foam -> Effective thermal properties ------------------//
    Info<< "Constructing field rhoEff\n" << endl;
    volScalarField rhoEff
    (
        IOobject
        (
            "rhoEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rhoEff", dimDensity, SMALL)
    );

    Info<< "Constructing field CpEff\n" << endl;
    volScalarField CpEff
    (
        IOobject
        (
            "CpEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("CpEff", dimSpecificHeatCapacity, SMALL)
    );

    Info<< "Constructing field kappaEff\n" << endl;
    volScalarField kappaEff
    (
        IOobject
        (
            "kappaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaEff", dimThermalConductivity, SMALL)
    );
//------------------- NF4Foam -> Radiation ------------------//
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );
