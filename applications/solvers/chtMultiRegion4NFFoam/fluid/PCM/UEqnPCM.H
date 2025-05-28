    Info << "    Type: PCM" << endl;

    const dimensionedVector& g = gFluid[i];
    dimensionedScalar& TRef = TRefFluid[i];
    volScalarField& rhok = rhokFluid[i];
    volScalarField& gamma = gammaFluid[i];
    dimensionedScalar& Cu = CuFluid[i];
    dimensionedScalar& b = bFluid[i];

    volScalarField hm = nanoFluidModelPtr.hm();
    volScalarField Tm = nanoFluidModelPtr.Tm();
    volScalarField Tmr = nanoFluidModelPtr.Tmr();


    // Solve the Momentum equation
    rhok = 1.0 - beta*(T - TRef);
    gamma = 0.5 * erf( 4.0 * ( T - Tm ) / Tmr ) + scalar(0.5);
    volScalarField SmPCM("SmPCM", Cu * pow( ( 1 - gamma ) , 2 ) / ( pow( gamma , 3 ) + b));
    Info << "max(p): " << max(p) << endl;


    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turb.divDevReff()
      + (1/rho)*fvm::SuSp(SmPCM, U)
    );

    UEqn().relax();

    pZones.addResistance(UEqn());

    if (momentumPredictor)
    {
         solve
         (
            UEqn()
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rhok)*(g & mesh.Sf())
              - fvc::snGrad(p)*mesh.magSf()
            )
        );
    }
