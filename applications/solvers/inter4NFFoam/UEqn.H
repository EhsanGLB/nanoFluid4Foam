
    surfaceScalarField alpha1f = min(max(fvc::interpolate(alpha1), scalar(0)), scalar(1));
    surfaceScalarField muf = alpha1f*fvc::interpolate(mu1) + (scalar(1) - alpha1f)*fvc::interpolate(mu2);

    surfaceScalarField muEff
    (
        "muEff",
        muf
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    pZones.addResistance(UEqn);
    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(pd)
                )*mesh.magSf()
            )
        );
    }
