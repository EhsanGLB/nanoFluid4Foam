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
#include "Krieger5.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace baseFluidDynamicViscosityModels
{
    defineTypeNameAndDebug(Krieger5, 0);
    addToRunTimeSelectionTable(baseFluidDynamicViscosityModel, Krieger5, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Krieger5::Krieger5
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    baseFluidDynamicViscosityModel(nanoFluidPropertiesDict, U, p, T),
    Krieger5Coeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    a_(Krieger5Coeffs_.lookup("a")),
    b_(Krieger5Coeffs_.lookup("b")),
    c_(Krieger5Coeffs_.lookup("c")),
    beta_(Krieger5Coeffs_.lookup("beta")),
    lambda_(Krieger5Coeffs_.lookup("lambda")),
    nuK_(Krieger5Coeffs_.lookup("nuK")),
    muPlasma_(Krieger5Coeffs_.lookup("muPlasma")),
    Hcrit_(Krieger5Coeffs_.lookup("Hcrit"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField Krieger5::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));
    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    mu_ = muPlasma_ * pow( (1-H/Hcrit_) , -1*((a_+b_*exp(-1*c_*H)+beta_*pow((1+pow(lambda_*strainRate(),2)),(-1*nuK_)))));

    return mu_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace baseFluidDynamicViscosityModels
} // End namespace Foam

// ************************************************************************* //
