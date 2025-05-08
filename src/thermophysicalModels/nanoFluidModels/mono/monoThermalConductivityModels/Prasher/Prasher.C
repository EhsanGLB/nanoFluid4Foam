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

#include "Prasher.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(Prasher, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, Prasher, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Prasher::Prasher
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
    PrasherCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    Aeco_(PrasherCoeffs_.lookup("Aeco")),
    Meco_(readScalar(PrasherCoeffs_.lookup("Meco"))),
    Rb_(PrasherCoeffs_.lookup("Rb")),
    dnp_(PrasherCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Prasher::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField rhonp_ = particleModelPtr_->rho();
    volScalarField kappanp_ = particleModelPtr_->kappa();


    volScalarField ReB_ = ( rhobf_ / mubf_) * sqrt ( ( 18.0 * kB_ * T_ ) / ( mathematicalConstant::pi * rhonp_ * dnp_ ) );
    volScalarField Prbf_ = (mubf_ * Cpbf_) / kappabf_;
    volScalarField kmtx_ = kappabf_ * ( 1 + 0.25 * ReB_ * Prbf_ );
    volScalarField Binp_ = ( 2.0 * Rb_ * kmtx_ ) / dnp_;
    kappa_ = ( ( kappanp_ * ( 1 + 2 * Binp_ ) + 2 * kmtx_ ) + 2 * alpha_ * ( kappanp_ * ( 1 - Binp_ ) - kmtx_ ) ) / ( ( kappanp_ * ( 1 + 2 * Binp_ ) + 2 * kmtx_ ) - alpha_ * ( kappanp_ * ( 1 - Binp_ ) - kmtx_ ) ) * ( 1 + Aeco_ * pow( ReB_ , Meco_) * pow( Prbf_ , 0.333) * alpha_) * kappabf_ ;

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
