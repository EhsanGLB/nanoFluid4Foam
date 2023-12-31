/*
            volScalarField rUA = 1.0/UEqn.A();

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf())
                + fvc::ddtPhiCorr(rUA, U, phi);

            adjustPhi(phi, U, p);

            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
*/
    p.boundaryField().updateCoeffs();

    volScalarField rAU = 1.0/UEqn().A();
    U = rAU*UEqn().H();
    UEqn.clear();

    // Calculate under-relaxation consistent flux
    simple.calcSteadyConsistentFlux(phi, U);
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                fvc::interpolate(rAU)/simple.aCoeff(U.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
         ==
            fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Reconstruct consistent velocity after pressure equation
    simple.reconstructSteadyVelocity(U, rAU, p);
