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
#include "Yeleswarapu.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{
    defineTypeNameAndDebug(Yeleswarapu, 0);
    addToRunTimeSelectionTable(baseFluidDynamicViscosityModel, Yeleswarapu, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Yeleswarapu::Yeleswarapu
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    baseFluidDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T),
    YeleswarapuCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    a1_(YeleswarapuCoeffs_.lookup("a1")),
    a2_(YeleswarapuCoeffs_.lookup("a2")),
    a3_(YeleswarapuCoeffs_.lookup("a3")),
    b1_(YeleswarapuCoeffs_.lookup("b1")),
    b2_(YeleswarapuCoeffs_.lookup("b2")),
    b3_(YeleswarapuCoeffs_.lookup("b3")),
    k_(YeleswarapuCoeffs_.lookup("k")),
    muPlasma_(YeleswarapuCoeffs_.lookup("muPlasma"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Yeleswarapu::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));

    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    mu_ = muPlasma_ * (1-H) + H * 
            muPlasma_/muPlasma_.value() *
            ((b3_ * pow(H,3) + b2_ * pow(H,2) + b1_ * H) + 
            ((a3_ * pow(H,3) + a2_ * pow(H,2) + a1_ * H) - 
            (b3_ * pow(H,3) + b2_ * pow(H,2) + b1_ * H)) * 
            (1+log(1+k_*strainRate()))/(1+k_*strainRate()));

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
