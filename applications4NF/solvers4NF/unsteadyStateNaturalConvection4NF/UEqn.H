    //-NFFoam
    nuEff = turbulence->nut() + mu/rho;
    nuEff.correctBoundaryConditions();

    fvVectorMatrix UEqn
    (
          fvm::ddt(U)
        + fvm::div(phi, U)
        - fvm::laplacian(nuEff, U)//-nanoFluid4Foam
        - fvc::div(nuEff*dev(fvc::grad(U)().T()))//-nanoFluid4Foam
          ==
        -beta*(T - TRef)*g
    );

    solve(UEqn == -fvc::grad(p));
