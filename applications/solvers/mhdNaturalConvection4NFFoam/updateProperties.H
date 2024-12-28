    rho = nanoFluidModelPtr->rho();
    kappa = nanoFluidModelPtr->kappa();
    Cp = nanoFluidModelPtr->Cp();
    mu = nanoFluidModelPtr->mu();
    beta = nanoFluidModelPtr->beta();
    rhoR = nanoFluidModelPtr->rhoR();
    nanoFluidModelPtr->print();

    forAll(moltensPtr, i)
    {
        moltensPtr[i] = particlesProperties[i].molten();
    }

    rhok = 1.0 - beta*(T - TRef);


    turbulence->correct();
    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();


    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;

    DB = rhoR/mu0;
    DBU = 1.0/(2.0*mu0*rho);

