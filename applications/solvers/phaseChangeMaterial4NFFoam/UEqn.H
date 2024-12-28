    // Solve the momentum equation
    gamma = 0.5*erf(4.0*(T-Tm)/Tmr)+scalar(0.5);
    volScalarField SmPCM("SmPCM", Cu*pow((1.0-gamma), 2.0)/(pow(gamma, 3.0)+b));

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
      + (1/rho)*fvm::SuSp(SmPCM, U)
    );

    UEqn.relax();

    pZones.addResistance(UEqn);

    if (piso.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(rhok)*(g & mesh.Sf())
                  - fvc::snGrad(p)*mesh.magSf()
                )
            )
        );
    }
