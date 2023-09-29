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

#include "JangChoi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(JangChoi, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, JangChoi, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

JangChoi::JangChoi
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
    JangChoiCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    beta_(readScalar(JangChoiCoeffs_.lookup("beta"))),
    C_(readScalar(JangChoiCoeffs_.lookup("C"))),
    lbf_(JangChoiCoeffs_.lookup("lbf")),
    dbf_(JangChoiCoeffs_.lookup("dbf")),
    dnp_(JangChoiCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField JangChoi::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField kappanp_ = particleModelPtr_->kappa();

    volScalarField CRM_ = (kB_ * T_) / (3 * mathematicalConstant::pi * mubf_ * dnp_ * lbf_) ;
    volScalarField Renp_ = ( rhobf_ * CRM_ * dnp_ ) / mubf_;
    volScalarField Prbf_ = (mubf_ * Cpbf_) / kappabf_;
    kappa_ = kappabf_ * (1 - alpha_) + beta_ * kappanp_ * alpha_ + C_ * kappabf_ * (dbf_ / dnp_) * alpha_ * Prbf_ *pow(Renp_, 2);

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
