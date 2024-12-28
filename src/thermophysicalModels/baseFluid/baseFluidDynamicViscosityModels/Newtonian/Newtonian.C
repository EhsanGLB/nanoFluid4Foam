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

#include "addToRunTimeSelectionTable.H"

#include "dimensionSets4NFFoam.H"
#include "Newtonian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{
    defineTypeNameAndDebug(Newtonian, 0);
    addToRunTimeSelectionTable(baseFluidDynamicViscosityModel, Newtonian, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Newtonian::Newtonian
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    baseFluidDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T),
    baseFluidComponents_(nanoFluidPropertiesDict.lookup("baseFluid")),
    moleFraction_(nanoFluidPropertiesDict.lookup("moleFraction")),
    materialMixturePtr_(materialMixture::New(baseFluidComponents_, nanoFluidPropertiesDict))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Newtonian::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        mu_.internalField()[i] = materialMixturePtr_->mu(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            mu_.boundaryField()[i][j] = materialMixturePtr_->mu(pij_, Tij_, moleFraction_);
	}
    }

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
