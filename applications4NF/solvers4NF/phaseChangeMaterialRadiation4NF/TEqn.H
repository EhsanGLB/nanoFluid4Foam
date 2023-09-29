        // Solve energy equation
    kappaEff = (turbulence->nut()/Prt) * rho * Cp + kappa;
    kappaEff.correctBoundaryConditions();
    volScalarField T_ = T;

    fvScalarMatrix TEqn
    (
          fvm::ddt(rho*Cp, T)
        + fvm::div(phi*fvc::interpolate(rho*Cp), T)
        - fvm::laplacian(kappaEff, T)
        + rho*hm*4.0*exp(-pow(4.0*(T-Tm)/Tmr,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr*fvm::ddt(T)
        + rho*hm*4.0*exp(-pow(4.0*(T-Tm)/Tmr,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr*(U & fvc::grad(T))
==
          radiation->Ru()
        - fvm::Sp(4.0*radiation->Rp()*pow(T_,3), T)
        - radiation->Rp()*pow(T_,3)*(T_ - 4.0*T)
    );

    TEqn.relax();
    TEqn.solve();
    //thermo.correct();
    radiation->correct();
    rhoBB = rho*(scalar(1) - beta*(T - TRef));
