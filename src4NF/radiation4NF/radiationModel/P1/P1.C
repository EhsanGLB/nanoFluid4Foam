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

#include "P1.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(P1, 0);
        addToRadiationRunTimeSelectionTables(P1);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::P1::P1(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            T.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    )
{}


Foam::radiation::P1::P1(const dictionary& dict, const volScalarField& T)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            T.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::P1::~P1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::P1::read()
{
    if (radiationModel::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::P1::calculate()
{
    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();
    const volScalarField sigmaEff(scatter_->sigmaEff());

    const dimensionedScalar a0 ("a0", a_.dimensions(), ROOTVSMALL);

    // Construct diffusion
    const volScalarField gamma
    (
        IOobject
        (
            "gammaRad",
            G_.mesh().time().timeName(),
            G_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0/(3.0*a_ + sigmaEff + a0)
    );

    // Solve G transport equation
    solve
    (
        fvm::laplacian(gamma, G_)
      - fvm::Sp(a_, G_)
     ==
      - 4.0*(e_*radiation::sigmaSB*pow4(T_) + E_)
    );

    // Calculate radiative heat flux on boundaries.
    forAll(mesh_.boundaryMesh(), patchI)
    {
        Qr_.boundaryField()[patchI] =
            -gamma.boundaryField()[patchI]*G_.boundaryField()[patchI].snGrad();
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiation::P1::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*absorptionEmission_->eCont()*radiation::sigmaSB
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::P1::Ru() const
{
    const DimensionedField<scalar, volMesh>& G =
        G_.dimensionedInternalField();
    const DimensionedField<scalar, volMesh> E =
        absorptionEmission_->ECont()().dimensionedInternalField();
    const DimensionedField<scalar, volMesh> a =
        absorptionEmission_->aCont()().dimensionedInternalField();

    return a*G - 4.0*E;
}


// ************************************************************************* //
