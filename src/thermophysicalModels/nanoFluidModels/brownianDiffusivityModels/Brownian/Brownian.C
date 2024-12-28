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

#include "dimensionSets4NFFoam.H"//-nanoFluid4Foam
#include "Brownian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace brownianDiffusivityModels
{
    defineTypeNameAndDebug(Brownian, 0);
    addToRunTimeSelectionTable(brownianDiffusivityModel, Brownian, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Brownian::Brownian
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& alphasPtr,
    const autoPtr<baseFluid>& baseFluidPtr,
    const PtrList<particleModel>& particlesProperties
)
:
    brownianDiffusivityModel(nanoFluidPropertiesDict, U, p, T, alphasPtr, baseFluidPtr, particlesProperties),
    BrownianCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    dnp_(BrownianCoeffs_.lookup("dnp")),
    mubf_(baseFluidPtr_->mu()),
    DBs_(particlesProperties.size())
{

    forAll(DBs_, i)
    {
        DBs_.set
        (
            i,
            new volScalarField (IOobject(("DB_"+particlesProperties_[i].particleName()), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("DB_"+particlesProperties_[i].particleName()), dimDiffusivity, SMALL))
        );
    }

    forAll(DBs_, i)
    {
        DBs_[i] = kB_ * T_ / ( 3 * mathematicalConstant::pi * mubf_ * dnp_[i]*dd_);
    }

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
/*
inline const PtrList<volScalarField>& Brownian::DBs() const
{
    return DBs_;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace brownianDiffusivityModels
} // End namespace Foam

// ************************************************************************* //
