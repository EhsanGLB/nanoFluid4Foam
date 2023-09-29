    //-nanoFluid4Foam
    nuEff = turbulence->nut() + mu/rho;
    nuEff.correctBoundaryConditions();

    tmp<fvVectorMatrix> UEqn
    (
          fvm::div(phi, U)
        - fvm::laplacian(nuEff, U)//-nanoFluid4Foam
        - fvc::div(nuEff*dev(fvc::grad(U)().T()))//-nanoFluid4Foam
    );

    UEqn().relax();

    // Momentum solution
    solve(relax(UEqn()) == -fvc::grad(p));
