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
#include "HerschelBulkley.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{
    defineTypeNameAndDebug(HerschelBulkley, 0);
    addToRunTimeSelectionTable(baseFluidDynamicViscosityModel, HerschelBulkley, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HerschelBulkley::HerschelBulkley
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    baseFluidDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T),
    HerschelBulkleyCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    k_(HerschelBulkleyCoeffs_.lookup("k")),
    n_(HerschelBulkleyCoeffs_.lookup("n")),
    tau0_(HerschelBulkleyCoeffs_.lookup("tau0")),
    mu0_(HerschelBulkleyCoeffs_.lookup("mu0"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField HerschelBulkley::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));

    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    dimensionedScalar rh("rh", dimDensity, 1.0);
    tmp<volScalarField> sr(strainRate());

    mu_ = (
min
(
mu0_,
(tau0_ + k_* rtone *( pow(tone * sr(), n_) - pow(tone*tau0_/mu0_, n_)))
/
(max(sr(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)))
)
);

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
