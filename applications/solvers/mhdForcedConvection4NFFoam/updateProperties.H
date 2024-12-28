    rho = nanoFluidModelPtr->rho();
    kappa = nanoFluidModelPtr->kappa();
    Cp = nanoFluidModelPtr->Cp();
    mu = nanoFluidModelPtr->mu();
    rhoR = nanoFluidModelPtr->rhoR();
    nanoFluidModelPtr->print();

    forAll(moltensPtr, i)
    {
        moltensPtr[i] = particlesProperties[i].molten();
    }


    turbulence->correct();
    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();


    rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;


    DB = rhoR/mu0;
    DBU = 1.0/(2.0*mu0*rho);







