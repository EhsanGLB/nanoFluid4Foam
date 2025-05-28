    const fvMesh& mesh = fluidRegions[i];

    volScalarField& T = TFluid[i];
    volScalarField& p = pFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];


    volScalarField& rho = rhoFluid[i];
    volScalarField& Cp = CpFluid[i];
    volScalarField& kappa = kappaFluid[i];
    volScalarField& mu = muFluid[i];
    volScalarField& beta = betaFluid[i];
    word& type = typeFluid[i];



    nanoFluidModel& nanoFluidModelPtr = nanoFluidModelPtrFluid[i];
    nanoFluidModel& nanoFluidModelObj = nanoFluidModelPtr;
    rho = nanoFluidModelPtr.rho();
    kappa = nanoFluidModelPtr.kappa();
    Cp = nanoFluidModelPtr.Cp();
    mu = nanoFluidModelPtr.mu();
    rho.correctBoundaryConditions();
    Cp.correctBoundaryConditions();
    kappa.correctBoundaryConditions();
    mu.correctBoundaryConditions();

    if(type == "naturalConvection")
    {
        beta = nanoFluidModelPtr.beta();
        beta.correctBoundaryConditions();
    }

    PtrList<volScalarField>& moltensPtr = moltensPtrFluid[i];
    const PtrList<particleModel>& particlesProperties = nanoFluidModelObj.particlesProperties();
    forAll(moltensPtr, j)
    {
        moltensPtr[j] = particlesProperties[j].molten();
    }


    mesh.schemesDict().setFluxRequired(p.name());




    incompressible::turbulenceModel& turb = turbulence[i];
    dimensionedScalar& Prt = PrtFluid[i];
    turb.correct();
    volScalarField& kappat = kappatFluid[i];
    kappat = (turb.nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();


    porousZones& pZones = pZonesFluid[i];
    volScalarField& porosity = porosityFluid[i];
    volScalarField& rhoPM = rhoPMFluid[i];
    volScalarField& CpPM = CpPMFluid[i];
    volScalarField& kappaPM = kappaPMFluid[i];
    volScalarField& hsPM = hsPMFluid[i];


    volScalarField& rhoEff = rhoEffFluid[i];
    volScalarField& CpEff = CpEffFluid[i];
    volScalarField& kappaEff = kappaEffFluid[i];


    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;
    rhoEff.correctBoundaryConditions();
    CpEff.correctBoundaryConditions();
    kappaEff.correctBoundaryConditions();

    radiation::radiationModel& rad = radiation[i];
