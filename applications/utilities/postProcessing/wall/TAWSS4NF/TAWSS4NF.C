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
    TAWSS4NF

Description
    Calculates and reports wall shear stress for all patches, for the
    specified times when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"

#include "basicPsiThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& WSS4NF
)
{
    IOobject muHeader
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!muHeader.headerOk())
    {
        Info<< "    no mu field" << endl;
        return;
    }

    Info<< "Reading field mu\n" << endl;
    volScalarField mu(muHeader, mesh);

    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> model
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    //const volSymmTensorField Reff(model->devReff());
    const volSymmTensorField Reff( -mu*dev(twoSymm(fvc::grad(U))) );

    forAll(WSS4NF.boundaryField(), patchI)
    {
        WSS4NF.boundaryField()[patchI] =
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI];
    }
}


void calcCompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& WSS4NF
)
{
    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.headerOk())
    {
        Info<< "    no rho field" << endl;
        return;
    }

    Info<< "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

    #include "compressibleCreatePhi.H"

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> model
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const volSymmTensorField Reff(model->devRhoReff());

    forAll(WSS4NF.boundaryField(), patchI)
    {
        WSS4NF.boundaryField()[patchI] =
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI];
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::validOptions.insert("compressible","");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    bool compressible = args.optionFound("compressible");

    volVectorField TAWSS4NF
    (
        IOobject
        (
            "TAWSS4NF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "TAWSS4NF",
            dimMass/dimLength/sqr(dimTime),
            vector::zero
        )
    );

    int nTime = 1;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volVectorField WSS4NF
        (
            IOobject
            (
                "WSS4NF",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "WSS4NF",
                dimMass/dimLength/sqr(dimTime),
                vector::zero
            )
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

            if (compressible)
            {
                calcCompressible(mesh, runTime, U, WSS4NF);
            }
            else
            {
                calcIncompressible(mesh, runTime, U, WSS4NF);
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        /*Info<< "Writing wall shear stress to field " << WSS4NF.name()
            << nl << endl;
        WSS4NF.write();*/



        TAWSS4NF += WSS4NF;

        TAWSS4NF /= nTime;
        Info<< "Writing average wall shear stress to field " << TAWSS4NF.name()
            << nl << endl;
        TAWSS4NF.write();
        TAWSS4NF *= nTime;
        nTime++;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
