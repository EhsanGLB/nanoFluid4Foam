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

#include "Xuan.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(Xuan, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, Xuan, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Xuan::Xuan
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
    XuanCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    lambda_(readScalar(XuanCoeffs_.lookup("lambda"))),
    Df_(readScalar(XuanCoeffs_.lookup("Df"))),
    C_(readScalar(XuanCoeffs_.lookup("C"))),
    Nnp_(readInt(XuanCoeffs_.lookup("Nnp"))),
    c_(XuanCoeffs_.lookup("c")),
    dnp_(XuanCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Xuan::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField rhonp_ = particleModelPtr_->rho();
    volScalarField kappanp_ = particleModelPtr_->kappa();
    volScalarField Cpnp_ = particleModelPtr_->Cp();

    dimensionedScalar rc_ = 0.5 * dnp_ * exp( ( scalar(Nnp_) - C_ ) / Df_ );
    kappa_ = ( ( kappanp_ + ( lambda_ - 1 ) * kappabf_ - ( lambda_ - 1 ) * ( kappabf_ - kappanp_ ) * alpha_ ) / ( kappanp_ + ( lambda_ - 1 ) * kappabf_ + ( kappabf_ - kappanp_ ) * alpha_ ) + ( rhonp_ * Cpnp_* alpha_ / ( 2 * kappabf_) ) * sqrt( (c_ * kB_ * T_) / (3 * mathematicalConstant::pi * mubf_ * rc_ ) ) ) * kappabf_;

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
