    //-NFFoam
    nuEff = turbulence->nut() + mu/rho;
    nuEff.correctBoundaryConditions();

    tmp<fvVectorMatrix> UEqn
    (
          fvm::div(phi, U)
        - fvm::laplacian(nuEff, U)//-nanoFluid4Foam
        - fvc::div(nuEff*dev(fvc::grad(U)().T()))//-nanoFluid4Foam
          ==
        - beta*(1.0 - alphasPtr[0])*(T - TRef)*g
        + ((rhop - rho)/rho)*betaC*(alphasPtr[0] - CRef)*g
        + ((rhom - rho)/rho)*betaM*(alpham - NRef)*g
    );

    solve(relax(UEqn()) == -fvc::grad(p));
