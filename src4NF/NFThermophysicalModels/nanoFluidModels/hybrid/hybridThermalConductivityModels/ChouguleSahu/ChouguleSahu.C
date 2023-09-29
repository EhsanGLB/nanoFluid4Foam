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

#include "ChouguleSahu.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace hybridThermalConductivityModels
{
    defineTypeNameAndDebug(ChouguleSahu, 0);
    addToRunTimeSelectionTable(hybridThermalConductivityModel, ChouguleSahu, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ChouguleSahu::ChouguleSahu
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
    hybridThermalConductivityModel(nanoFluidPropertiesDict, U, p, T, alphasPtr, baseFluidPtr, particlesProperties),
    ChouguleSahuCoeffs_(nanoFluidPropertiesDict.subDict(typeName + "Coeffs")),
    dbf_(ChouguleSahuCoeffs_.lookup("dbf")),
    dnpList_(ChouguleSahuCoeffs_.lookup("dnp"))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField ChouguleSahu::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField alpha_(IOobject("alpha", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("alpha", dimless, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();

    forAll(particlesComponents_, i)
    {
        alpha_ += alphasPtr_[i];
    }

    kappa_ = kappabf_;
    forAll(particlesComponents_, i)
    {
        kappa_ += ( ( alphasPtr_[i] * particlesProperties_[i].kappa() * dbf_.value() ) / ( ( 1 - alpha_ ) * dnpList_[i] ) );
    }

    return kappa_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hybridThermalConductivityModel
} // End namespace Foam

// ************************************************************************* //
