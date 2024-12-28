{
	tmp<fv::convectionScheme<scalar> > mvConvection
	(
	    fv::convectionScheme<scalar>::New
	    (
		mesh,
		fields2,
		phi,
		mesh.schemesDict().divScheme("div(phi,alphai)")
	    )
	);



	forAll(alphasPtr2, i)
	{
	    volScalarField& alphai = alphasPtr2[i];
	    volScalarField DBi("DB_" + particlesComponents2[i][0], DBs2[i]);
	    volScalarField DBEffi ("DBEff_" + particlesComponents2[i][0], turbulence->nut()/Sct + DBi);
	    DBEffi.correctBoundaryConditions();
	    volScalarField DTi("DT_" + particlesComponents2[i][0], DTs2[i]);
	    solve
	    (
		  fvm::ddt(alphai)
		+ mvConvection->fvmDiv(phi, alphai)
		- fvm::laplacian(DBEffi, alphai)
		- fvc::laplacian((DTi/T), T)
	    );
	}


	forAll(alphasPtr2, i)
	{
	    cout << boldBlue << "The minimum and maximum of " << alphasPtr2[i].name() << ": " << boldWhite << min(gMin(alphasPtr2[i].internalField()), gMin(alphasPtr2[i].boundaryField())) << " ," << max(gMax(alphasPtr2[i].internalField()), gMax(alphasPtr2[i].boundaryField())) << reset << "\n";
	}
}
