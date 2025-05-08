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

#include "WalvekarCNT.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(WalvekarCNT, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, WalvekarCNT, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WalvekarCNT::WalvekarCNT
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
    WalvekarCNTCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    lcnt_(WalvekarCNTCoeffs_.lookup("lcnt")),
    T0_(WalvekarCNTCoeffs_.lookup("T0")),
    dcnt_(WalvekarCNTCoeffs_.lookup("dcnt")),
    dbf_(WalvekarCNTCoeffs_.lookup("dbf"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField WalvekarCNT::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField kappanp_ = particleModelPtr_->kappa();

    dimensionedScalar C_ = ( 85.0 * pow( kB_ , 2.0 ) ) / ( 72.0 * pow( mathematicalConstant::pi , 2.0) );
    kappa_ = kappabf_ * ( 1.0 + ( 2.0 * kappanp_ * dbf_ * alpha_ * ( 0.5 * dcnt_ + lcnt_ ) ) / ( 3.0 * dcnt_ * lcnt_ * kappabf_ * ( 1.0 - alpha_ ) ) ) + log( lcnt_ / dcnt_ ) * ( C_ * alpha_ * ( T_ - T0_ ) ) / ( pow( 0.5 * dcnt_ , 2.0) * pow( lcnt_ , 2.0) * mubf_ );
    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
