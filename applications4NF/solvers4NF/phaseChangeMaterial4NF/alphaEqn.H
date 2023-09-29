tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemesDict().divScheme("div(phi,alphai)")
    )
);



forAll(alphasPtr, i)
{
    volScalarField& alphai = alphasPtr[i];
    volScalarField DBi("DB_" + particlesComponents[i][0], DBs[i]);
    volScalarField DBEffi ("DBEff_" + particlesComponents[i][0], turbulence->nut()/Sct + DBi);
    DBEffi.correctBoundaryConditions();
    volScalarField DTi("DT_" + particlesComponents[i][0], DTs[i]);
    solve
    (
        mvConvection->fvmDiv( phi , alphai )
        - fvm::laplacian( DBEffi , alphai )
        - fvc::laplacian( (DTi/T) , T )
    );
}


forAll(alphasPtr, i)
{
    cout << boldBlue << "The minimum and maximum of " << alphasPtr[i].name() << ": " << boldWhite << min(gMin(alphasPtr[i].internalField()), gMin(alphasPtr[i].boundaryField())) << " ," << max(gMax(alphasPtr[i].internalField()), gMax(alphasPtr[i].boundaryField())) << reset << "\n";
}
