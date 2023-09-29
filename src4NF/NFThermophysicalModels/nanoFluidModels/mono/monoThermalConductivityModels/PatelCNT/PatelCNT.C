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

#include "PatelCNT.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(PatelCNT, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, PatelCNT, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PatelCNT::PatelCNT
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
    monoThermalConductivityModel(nanoFluidPropertiesDict, U, p, T, alpha, baseFluidPtr, particleModelPtr),
    PatelCNTCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    dbf_(PatelCNTCoeffs_.lookup("dbf")),
    dcnt_(PatelCNTCoeffs_.lookup("dcnt"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField PatelCNT::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField kappanp_ = particleModelPtr_->kappa();

    kappa_ = ( 1.0 + ( kappanp_ * alpha_ * dbf_ ) / ( kappabf_ * ( 1.0 - alpha_ ) * dcnt_ ) ) * kappabf_;

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
