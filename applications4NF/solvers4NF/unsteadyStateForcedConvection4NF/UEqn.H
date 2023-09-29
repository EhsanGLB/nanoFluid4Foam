//-nanoFluid4Foam
    nuEff = turbulence->nut() + mu/rho;
    nuEff.correctBoundaryConditions();

// Time derivative matrix
fvVectorMatrix ddtUEqn(fvm::ddt(U));

// Convection-diffusion matrix
fvVectorMatrix HUEqn
(
      fvm::div(phi, U)
    - fvm::laplacian(nuEff, U)//-nanoFluid4Foam
    - fvc::div(nuEff*dev(fvc::grad(U)().T()))//-nanoFluid4Foam
);

if (pimple.momentumPredictor())
{
    // Solve momentum predictor
    solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
}
else
{
    // Explicit update
    U = (ddtUEqn.H() + HUEqn.H() - fvc::grad(p))/(HUEqn.A() + ddtUEqn.A());
    U.correctBoundaryConditions();
}
