{
    surfaceScalarField alpha1f = min(max(fvc::interpolate(alpha1), scalar(0)), scalar(1));
    surfaceScalarField kappaf = alpha1f*fvc::interpolate(kappa1) + (1.0 - alpha1f)*fvc::interpolate(kappa2);

    //rhoCp = alpha1*rho1*Cp1 + (scalar(1) - alpha1)*rho2*Cp2;
    //kappa = alpha1*kappa1 + (scalar(1) - alpha1)*kappa2;
    volScalarField T_ = T;


    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaf, T)
      ==
        radiation->Ru()
      - fvm::Sp(4.0*radiation->Rp()*pow(T_,3), T)
      - radiation->Rp()*pow(T_,3)*(T_ - 4.0*T)
    );

    TEqn.solve();

    radiation->correct();
}
