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

#include "Yang.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoThermalConductivityModels
{
    defineTypeNameAndDebug(Yang, 0);
    addToRunTimeSelectionTable(monoThermalConductivityModel, Yang, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Yang::Yang
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
    YangCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    Rb_(YangCoeffs_.lookup("Rb")),
    dnp_(YangCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Yang::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField rhonp_ = particleModelPtr_->rho();

    dimensionedScalar Vnp_ = ( mathematicalConstant::pi / 6.0) * pow( dnp_ , 3.0);
    volScalarField resis_ = dnp_ / ( 2 * Rb_ * kappabf_);
    volScalarField taunp_ = ( rhonp_ * Vnp_ ) / ( 3 * mathematicalConstant::pi * mubf_ * dnp_ );
    volScalarField unp_ = sqrt( ( 3 * kB_ * T_ ) / ( rhonp_ * Vnp_ ) );
    kappa_ = (1 + 3 * alpha_ * ( ( resis_ - 1 ) / ( resis_ + 2 ) ) ) * kappabf_ + 157.5 * alpha_ * taunp_ * pow( unp_ , 2.0 ) * ( rhobf_ * Cpbf_ );

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoThermalConductivityModels
} // End namespace Foam

// ************************************************************************* //
