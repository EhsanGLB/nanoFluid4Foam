    //------------ phase1 ------------//
    rho1 = nanoFluidModelPtr1->rho();
    kappa1 = nanoFluidModelPtr1->kappa();
    Cp1 = nanoFluidModelPtr1->Cp();
    mu1 = nanoFluidModelPtr1->mu();
    //nanoFluidModelPtr1->print();

    forAll(moltensPtr1, i)
    {
        moltensPtr1[i] = particlesProperties1[i].molten();
    }

    //------------ phase2 ------------//
    rho2 = nanoFluidModelPtr2->rho();
    kappa2 = nanoFluidModelPtr2->kappa();
    Cp2 = nanoFluidModelPtr2->Cp();
    mu2 = nanoFluidModelPtr2->mu();
    //nanoFluidModelPtr2->print();

    forAll(moltensPtr2, i)
    {
        moltensPtr2[i] = particlesProperties2[i].molten();
    }

    //------------ equivalent ------------//

    turbulence->correct();
    kappat = (turbulence->nut()/Prt)*rhoCp;
    kappat.correctBoundaryConditions();


    /*rhoEff = porosity*rho + (1-porosity)*rhoPM;
    CpEff = (porosity*rho*Cp + (1-porosity)*rhoPM*CpPM)/rhoEff;
    kappaEff = porosity*(kappat + kappa) + (1-porosity)*kappaPM;*/

