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

#include "helicity4NF.H"
#include "volFields.H"
#include "dictionary.H"
#include "fvcCurl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helicity4NF, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helicity4NF::helicity4NF
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    UName_("U"),
    outputName_(typeName)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "helicity4NF::helicity4NF"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* helicity4NFPtr
        (
            new volScalarField
            (
                IOobject
                (
                    outputName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimLength/sqr(dimTime), 0.0)
            )
        );

        mesh.objectRegistry::store(helicity4NFPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helicity4NF::~helicity4NF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helicity4NF::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        if (UName_ != "U")
        {
            outputName_ = typeName + "(" + UName_ + ")";
        }
    }
}


void Foam::helicity4NF::execute()
{
    if (active_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        volScalarField& helicity4NF =
            const_cast<volScalarField&>
            (
                obr_.lookupObject<volScalarField>(outputName_)
            );

        helicity4NF = U & (fvc::curl(U));
    }
}


void Foam::helicity4NF::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::helicity4NF::timeSet()
{
    // Do nothing
}


void Foam::helicity4NF::write()
{
    if (active_)
    {
        const volScalarField& helicity4NF =
            obr_.lookupObject<volScalarField>(outputName_);

        Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << helicity4NF.name() << nl
            << endl;

        helicity4NF.write();
    }
}


// ************************************************************************* //
