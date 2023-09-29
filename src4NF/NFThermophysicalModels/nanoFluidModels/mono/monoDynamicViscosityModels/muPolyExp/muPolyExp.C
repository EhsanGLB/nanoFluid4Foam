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

#include "NFDimensionSet.H"
#include "muPolyExp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace monoDynamicViscosityModels
{
    defineTypeNameAndDebug(muPolyExp, 0);
    addToRunTimeSelectionTable(monoDynamicViscosityModel, muPolyExp, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muPolyExp::muPolyExp
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
    muPolyExpCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    a0_(readScalar(muPolyExpCoeffs_.lookup("a0"))),
    a1_(readScalar(muPolyExpCoeffs_.lookup("a1"))),
    a2_(readScalar(muPolyExpCoeffs_.lookup("a2"))),
    a3_(readScalar(muPolyExpCoeffs_.lookup("a3"))),
    a4_(readScalar(muPolyExpCoeffs_.lookup("a4"))),
    a5_(readScalar(muPolyExpCoeffs_.lookup("a5"))),
    a6_(readScalar(muPolyExpCoeffs_.lookup("a6"))),
    b_(readScalar(muPolyExpCoeffs_.lookup("b")))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField muPolyExp::mu() const
{
    volScalarField mu_(IOobject("muNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muNF", dimDynamicViscosity, SMALL));
    volScalarField mubf_ = baseFluidPtr_->mu();
    dimensionedScalar dimlessT_("dimlessT", dimensionSet(0, 0, 0, 1, 0, 0, 0), scalar(1.0));

    mu_ = mubf_ * ( a6_ * pow( alpha_ , 6.0) + a5_ * pow( alpha_ , 5.0) + a4_ * pow( alpha_ , 4.0) + a3_ * pow( alpha_ , 3.0) + a2_ * pow( alpha_ , 2.0) + a1_ * alpha_ + a0_ ) * exp( b_ * (T_/dimlessT_) );

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace monoDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
