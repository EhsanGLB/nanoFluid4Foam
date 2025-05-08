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

#include "dimensionSets4NFFoam.H"
#include "Masoumi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoDynamicViscosityModels
{
    defineTypeNameAndDebug(Masoumi, 0);
    addToRunTimeSelectionTable(monoDynamicViscosityModel, Masoumi, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Masoumi::Masoumi
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
    MasoumiCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    c1_(MasoumiCoeffs_.lookup("c1")),
    c2_(MasoumiCoeffs_.lookup("c2")),
    c3_(MasoumiCoeffs_.lookup("c3")),
    c4_(MasoumiCoeffs_.lookup("c4")),
    dnp_(MasoumiCoeffs_.lookup("dnp"))
{
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Masoumi::mu() const
{
    volScalarField mu_(IOobject("muNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muNF", dimDynamicViscosity, SMALL));
    volScalarField mubf_ = baseFluidPtr_->mu();
    volScalarField rhonp_ = particleModelPtr_->rho();

    volScalarField VB_ = ( 1.0 / dnp_ ) * sqrt( ( 18.0 * kB_ * T_ ) / ( mathematicalConstant::pi * rhonp_ * dnp_ ) );
    volScalarField delta_ = dnp_ * pow( ( mathematicalConstant::pi / ( 6.0 * alpha_ ) ) , 0.333 );
    volScalarField C_ = ( 1.0 / mubf_ ) * ( ( c1_ * dnp_ + c2_) * alpha_ + ( c3_ * dnp_ + c4_ ) );
    mu_ = mubf_ * ( 1.0 + ( rhonp_ * VB_ * pow( dnp_ , 2.0) ) / ( 72.0 * C_ * delta_ * mubf_ ) );

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
