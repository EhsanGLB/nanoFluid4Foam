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

\*---------------------------------------------------------------------------*/

#include "magWSS4NFFoam.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "wallPolyPatch.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(magWSS4NFFoam, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::magWSS4NFFoam::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall shear stress");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    file() << endl;
}


void Foam::magWSS4NFFoam::calcWSS4NFFoam
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff,
    volScalarField& WSS4NFFoam
)
{
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        scalarField& ssp = WSS4NFFoam.boundaryField()[patchI];
        const vectorField& Sfp = mesh.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchI];
        const symmTensorField& Reffp = Reff.boundaryField()[patchI];

        ssp = mag( (-Sfp/magSfp) & Reffp );

        scalar minSsp = gMin(ssp);
        scalar maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            file() << mesh.time().value()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        if (log_) Info<< "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::magWSS4NFFoam::magWSS4NFFoam
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    patchSet_(),
    UName_("U"),
    muName_("mu"),
    rhoName_("rho")
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "magWSS4NFFoam::magWSS4NFFoam"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* magWSS4NFFoamPtr
        (
            new volScalarField
            (
                IOobject
                (
                    type(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "0",
                    dimMass/dimLength/sqr(dimTime),
                    0.0
                )
            )
        );

        mesh.objectRegistry::store(magWSS4NFFoamPtr);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::magWSS4NFFoam::~magWSS4NFFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::magWSS4NFFoam::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        muName_ = dict.lookupOrDefault<word>("muName", "mu");
        rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");

        log_ = dict.lookupOrDefault<Switch>("log", true);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        patchSet_ =
            mesh.boundaryMesh().patchSet
            (
                wordReList(dict.lookupOrDefault("patches", wordReList()))
            );

        Info<< type() << " " << name_ << ":" << nl;

        if (patchSet_.empty())
        {
            forAll(pbm, patchI)
            {
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    patchSet_.insert(patchI);
                }
            }

            Info<< "    processing all wall patches" << nl << endl;
        }
        else
        {
            Info<< "    processing wall patches: " << nl;
            labelHashSet filteredPatchSet;
            forAllConstIter(labelHashSet, patchSet_, iter)
            {
                label patchI = iter.key();
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    filteredPatchSet.insert(patchI);
                    Info<< "        " << pbm[patchI].name() << endl;
                }
                else
                {
                    WarningIn("void magWSS4NFFoam::read(const dictionary&)")
                        << "Requested wall shear stress on non-wall boundary "
                        << "type patch: " << pbm[patchI].name() << endl;
                }
            }

            Info<< endl;

            patchSet_ = filteredPatchSet;
        }
    }
}


void Foam::magWSS4NFFoam::execute()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    if (active_)
    {
        ////////
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField& mu = obr_.lookupObject<volScalarField>(muName_);
        const volScalarField& rho = obr_.lookupObject<volScalarField>(rhoName_);

        functionObjectFile::write();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField& magWSS4NFFoam =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(type())
            );

        if (log_) Info<< type() << " " << name_ << " output:" << nl;


        tmp<volSymmTensorField> Reff;
        tmp<volScalarField> nut;
        /*if (mesh.foundObject<cmpModel>("turbulenceProperties"))
        {
            const cmpModel& model =
                mesh.lookupObject<cmpModel>("turbulenceProperties");

            Reff = model.devRhoReff();
            //nut = model.mut();
        }
        else if (mesh.foundObject<icoModel>("turbulenceProperties"))
        {
            const icoModel& model =
                mesh.lookupObject<icoModel>("turbulenceProperties");

            Reff = model.devReff();
            //nut = model.nut();
            //Reff = -nut()*dev(twoSymm(fvc::grad(U)));
        }
        else
        {
            FatalErrorIn("void Foam::magWSS4NFFoam::execute()")
                << "Unable to find turbulence model in the "
                << "database" << exit(FatalError);
        }*/
            Reff = -mu*dev(twoSymm(fvc::grad(U)));

        calcWSS4NFFoam(mesh, Reff(), magWSS4NFFoam);
    }
}


void Foam::magWSS4NFFoam::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::magWSS4NFFoam::timeSet()
{
    // Do nothing
}


void Foam::magWSS4NFFoam::write()
{
    if (active_)
    {
        functionObjectFile::write();

        const volScalarField& magWSS4NFFoam =
            obr_.lookupObject<volScalarField>(type());

        if (log_) Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << magWSS4NFFoam.name() << nl
            << endl;

        magWSS4NFFoam.write();
    }
}


// ************************************************************************* //
