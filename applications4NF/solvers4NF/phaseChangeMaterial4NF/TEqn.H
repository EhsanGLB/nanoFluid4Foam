        // Solve energy equation
    kappaEff = (turbulence->nut()/Prt) * rho * Cp + kappa;
    kappaEff.correctBoundaryConditions();
    //fvScalarMatrix Sh = (1/Cp)*radiation->Sh(thermo);

    fvScalarMatrix TEqn
    (
            fvm::ddt(rho*Cp, T)
          + fvm::div(phi*fvc::interpolate(rho*Cp), T)
          - fvm::laplacian(kappaEff, T)
          + rho*hm*4.0*exp(-pow(4.0*(T-Tm)/Tmr,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr*fvm::ddt(T)
          + rho*hm*4.0*exp(-pow(4.0*(T-Tm)/Tmr,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr*(U & fvc::grad(T))
/*==
           (1/Cp)*radiation->Sh(thermo)*/
    );

    TEqn.relax();
    TEqn.solve();
    //thermo.correct();
    //radiation->correct();
    rhoBB = rho*(scalar(1) - beta*(T - TRef));
