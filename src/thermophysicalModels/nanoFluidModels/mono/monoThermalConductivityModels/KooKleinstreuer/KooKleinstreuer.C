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

#include "KooKleinstreuer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(KooKleinstreuer, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, KooKleinstreuer, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KooKleinstreuer::KooKleinstreuer
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
    KooKleinstreuerCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    lambda_(readScalar(KooKleinstreuerCoeffs_.lookup("lambda"))),
    Abeta_(KooKleinstreuerCoeffs_.lookup("Abeta")),
    Bbeta_(KooKleinstreuerCoeffs_.lookup("Bbeta")),
    Af_(KooKleinstreuerCoeffs_.lookup("Af")),
    Bf_(KooKleinstreuerCoeffs_.lookup("Bf")),
    Cf_(KooKleinstreuerCoeffs_.lookup("Cf")),
    Df_(KooKleinstreuerCoeffs_.lookup("Df")),
    dnp_(KooKleinstreuerCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField KooKleinstreuer::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField rhonp_ = particleModelPtr_->rho();
    volScalarField kappanp_ = particleModelPtr_->kappa();
    volScalarField Cpnp_ = particleModelPtr_->Cp();

    volScalarField beta_ = Abeta_ * pow( ( 100* alpha_ ) , Bbeta_);
    volScalarField f_ = ( Af_ * alpha_ + Bf_ ) * T_ + Cf_ * alpha_ + Df_;
    kappa_ = ( ( kappanp_ + ( lambda_ - 1 ) * kappabf_ - ( lambda_ - 1 ) * ( kappabf_ - kappanp_ ) * alpha_ ) / ( kappanp_ + ( lambda_ - 1 ) * kappabf_ + ( kappabf_ - kappanp_ ) * alpha_ ) + ( 5.0e4 * beta_ * rhonp_ * Cpnp_* alpha_ ) * sqrt( (kB_ * T_) / (rhonp_ * dnp_) ) * f_ ) * kappabf_;

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
