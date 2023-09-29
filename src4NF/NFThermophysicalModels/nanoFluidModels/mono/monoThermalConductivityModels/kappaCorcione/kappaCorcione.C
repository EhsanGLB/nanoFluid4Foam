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

#include "kappaCorcione.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(kappaCorcione, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, kappaCorcione, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kappaCorcione::kappaCorcione
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
    kappaCorcioneCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    Tfr_(kappaCorcioneCoeffs_.lookup("T_fr")),
    dnp_(kappaCorcioneCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField kappaCorcione::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField kappanp_ = particleModelPtr_->kappa();

    volScalarField Prbf_ = (mubf_ * Cpbf_) / kappabf_;
    volScalarField Renp_ = ( 2.0 * rhobf_ * kB_ * T_ ) / ( mathematicalConstant::pi * pow(mubf_ , 2) * dnp_);
    kappa_ = ( 1 + 4.4 * pow( alpha_ , 0.66 ) * pow( ( kappanp_ / kappabf_) , 0.03 )  * pow( ( T_ / Tfr_) , 10.0 ) * pow( Prbf_ , 0.66 ) * pow( Renp_ , 0.4 ) ) * kappabf_;

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
