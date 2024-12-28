{
	tmp<fv::convectionScheme<scalar> > mvConvection
	(
	    fv::convectionScheme<scalar>::New
	    (
		mesh,
		fields1,
		phi,
		mesh.schemesDict().divScheme("div(phi,alphai)")
	    )
	);



	forAll(alphasPtr1, i)
	{
	    volScalarField& alphai = alphasPtr1[i];
	    volScalarField DBi("DB_" + particlesComponents1[i][0], DBs1[i]);
	    volScalarField DBEffi ("DBEff_" + particlesComponents1[i][0], turbulence->nut()/Sct + DBi);
	    DBEffi.correctBoundaryConditions();
	    volScalarField DTi("DT_" + particlesComponents1[i][0], DTs1[i]);
	    solve
	    (
		  fvm::ddt(alphai)
		+ mvConvection->fvmDiv(phi, alphai)
		- fvm::laplacian(DBEffi, alphai)
		- fvc::laplacian((DTi/T), T)
	    );
	}


	forAll(alphasPtr1, i)
	{
	    cout << boldBlue << "The minimum and maximum of " << alphasPtr1[i].name() << ": " << boldWhite << min(gMin(alphasPtr1[i].internalField()), gMin(alphasPtr1[i].boundaryField())) << " ," << max(gMax(alphasPtr1[i].internalField()), gMax(alphasPtr1[i].boundaryField())) << reset << "\n";
	}
}
