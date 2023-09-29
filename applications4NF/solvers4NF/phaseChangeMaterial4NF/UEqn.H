    //-NFFoam
    nuEff = turbulence->nut() + mu/rho;
    nuEff.correctBoundaryConditions();
    volScalarField Vgamma ("Vgamma", Cu * pow( ( 1 - gamma ) , 2 ) / ( pow( gamma , 3 ) + b));

    fvVectorMatrix UEqn
    (
          fvm::ddt(U)
        + fvm::div(phi, U)
        - fvm::laplacian(nuEff, U)//-nanoFluid4Foam
        - fvc::div(nuEff*dev(fvc::grad(U)().T()))//-nanoFluid4Foam
        + (1/rho)*fvm::SuSp(Vgamma, U)//-nanoFluid4Foam
          ==
        -beta*(T - TRef)*g
    );

    solve(UEqn == -fvc::grad(p));
