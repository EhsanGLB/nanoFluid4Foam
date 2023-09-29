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

#include "NFDimensionSet.H"
#include "Mooney.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoDynamicViscosityModels
{
    defineTypeNameAndDebug(Mooney, 0);
    addToRunTimeSelectionTable(monoDynamicViscosityModel, Mooney, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Mooney::Mooney
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T,
    const volScalarField& alpha,
    const autoPtr<baseFluid>& baseFluidPtr,
    const autoPtr<particleModel>& particleModelPtr
)
:
    monoDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T, alpha, baseFluidPtr, particleModelPtr),
    MooneyCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    kesi_(readScalar(MooneyCoeffs_.lookup("kesi"))),
    zeta_(readScalar(MooneyCoeffs_.lookup("zeta")))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Mooney::mu() const
{
    volScalarField mu_(IOobject("muNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muNF", dimDynamicViscosity, SMALL));
    volScalarField mubf_ = baseFluidPtr_->mu();

    mu_ = mubf_ * exp( ( kesi_ * alpha_ ) / ( 1 - zeta_ * alpha_ ) );

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
