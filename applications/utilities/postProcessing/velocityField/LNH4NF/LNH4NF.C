/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    LNH4NF

Description
    Calculates and writes the LNH4NF of velocity field U.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);

        Info<< "    Calculating LNH4NF" << endl;

        volScalarField LNH4NF
        (
            IOobject
            (
                "LNH4NF",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
        );

        dimensionedScalar epsLNH("epsLNH",dimLength/sqr(dimTime), 1e-64);
        LNH4NF = ( U & (fvc::curl(U)) ) / (mag(U) * mag(fvc::curl(U)) + epsLNH );


        Info<< "LNH4NF max/min : "
            << max(LNH4NF).value() << " "
            << min(LNH4NF).value() << endl;

        if (writeResults)
        {
            LNH4NF.write();
            LNH4NF.write();
        }
    }
    else
    {
        Info<< "    No U" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
