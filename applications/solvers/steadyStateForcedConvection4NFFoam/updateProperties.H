    rho = nanoFluidModelPtr->rho();
    kappa = nanoFluidModelPtr->kappa();
    Cp = nanoFluidModelPtr->Cp();
    mu = nanoFluidModelPtr->mu();
    nanoFluidModelPtr->print();

    forAll(moltensPtr, i)
    {
        moltensPtr[i] = particlesProperties[i].molten();
    }


    turbulence->correct();
    kappat = (turbulence->nut()/Prt)*rho*Cp;
    kappat.correctBoundaryConditions();

    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;

