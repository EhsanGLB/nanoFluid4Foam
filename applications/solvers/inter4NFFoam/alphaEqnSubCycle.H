label nAlphaCorr
(
    readLabel(pimple.dict().lookup("nAlphaCorr"))
);

label nAlphaSubCycles
(
    readLabel(pimple.dict().lookup("nAlphaSubCycles"))
);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum = 0.0*rhoPhi;

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
#       include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
#   include "alphaEqn.H"
}

interface.correct();

rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
rhoCp = alpha1*rho1*Cp1 + (scalar(1) - alpha1)*rho2*Cp2;
