//-nanoFluid4Foam
{
    kappaEff = (turbulence->nut()/Prt) * rho * Cp + kappa;
    kappaEff.correctBoundaryConditions();

    fvScalarMatrix TEqn
    (
        fvm::div(phi*fvc::interpolate(rho*Cp), T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();
    TEqn.solve();
}
