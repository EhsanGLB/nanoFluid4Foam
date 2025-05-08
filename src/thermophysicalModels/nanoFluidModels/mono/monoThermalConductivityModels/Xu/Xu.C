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

#include "Xu.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(Xu, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, Xu, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Xu::Xu
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
    XuCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    lambda_(readScalar(XuCoeffs_.lookup("lambda"))),
    dd_(readScalar(XuCoeffs_.lookup("dd"))),
    Hco_(readScalar(XuCoeffs_.lookup("Hco"))),
    dnpmax_(XuCoeffs_.lookup("dnpmax")),
    dnpmin_(XuCoeffs_.lookup("dnpmin")),
    dnpavg_(XuCoeffs_.lookup("dnpavg")),
    Nunp_(XuCoeffs_.lookup("Nunp")),
    dbf_(XuCoeffs_.lookup("dbf"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Xu::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField kappanp_ = particleModelPtr_->kappa();

    volScalarField Prbf_ = (mubf_ * Cpbf_) / kappabf_;
    volScalarField Df_ = dd_ - log(alpha_) / log(dnpmax_ / dnpmin_);
    kappa_ = ( ( kappanp_ + ( lambda_ - 1 ) * kappabf_ - ( lambda_ - 1 ) * ( kappabf_ - kappanp_ ) * alpha_ ) / ( kappanp_ + ( lambda_ - 1 ) * kappabf_ + ( kappabf_ - kappanp_ ) * alpha_ ) + Hco_ * ( Nunp_ * dbf_ / Prbf_ ) * ( ( 2 - Df_ ) * Df_ / pow( (1 - Df_), 2.0) ) * ( pow( pow( ( dnpmax_ / dnpmin_ ) , ( 1 - Df_ ) ) - 1 , 2.0 ) / (pow( ( dnpmax_ / dnpmin_ ) , ( 2 - Df_ ) ) - 1 ) ) * ( 1 / dnpavg_) ) * kappabf_;

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
