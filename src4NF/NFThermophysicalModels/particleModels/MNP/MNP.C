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

#include "NFDimensionSet.H"//-nanoFluid4Foam
#include "MNP.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleModels
{
    defineTypeNameAndDebug(MNP, 0);
    addToRunTimeSelectionTable(particleModel, MNP, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MNP::MNP
(
    const wordList& particleList,
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    particleModel(particleList, nanoFluidPropertiesDict, U, p, T),
    MNPComponents_(particleSubDict_.lookup("Components")),
    MNPMoleFraction_(particleSubDict_.lookup("moleFraction")),
    materialMixturePtr_(materialMixture::New(MNPComponents_, nanoFluidPropertiesDict))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
const volScalarField MNP::rho() const
{
    volScalarField rho_(IOobject(("rho_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rho_"+particleName_), dimDensity, SMALL));

    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rho_.internalField()[i] = materialMixturePtr_->rho(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rho_.boundaryField()[i][j] = materialMixturePtr_->rho(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return rho_;
}


const volScalarField MNP::kappa() const
{
    volScalarField kappa_(IOobject(("kappa_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("kappa_"+particleName_), dimThermalConductivity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        kappa_.internalField()[i] = materialMixturePtr_->K(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            kappa_.boundaryField()[i][j] = materialMixturePtr_->K(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return kappa_;
}


const volScalarField MNP::Cp() const
{
    volScalarField Cp_(IOobject(("Cp_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Cp_"+particleName_), dimSpecificHeatCapacity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Cp_.internalField()[i] = materialMixturePtr_->cp(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Cp_.boundaryField()[i][j] = materialMixturePtr_->cp(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return Cp_;
}


const volScalarField MNP::beta() const
{
    volScalarField beta_(IOobject(("beta_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("beta_"+particleName_), dimThermalExpansion, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        beta_.internalField()[i] = materialMixturePtr_->beta(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            beta_.boundaryField()[i][j] = materialMixturePtr_->beta(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return beta_;
}


const volScalarField MNP::rhoR() const
{
    volScalarField rhoR_(IOobject(("rhoR_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("rhoR_"+particleName_), dimElectricalResistivity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rhoR_.internalField()[i] = materialMixturePtr_->rhoR(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rhoR_.boundaryField()[i][j] = materialMixturePtr_->rhoR(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return rhoR_;
}


const volScalarField MNP::mu() const
{
    volScalarField mu_(IOobject(("mu_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("mu_"+particleName_), dimDynamicViscosity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        mu_.internalField()[i] = materialMixturePtr_->mu(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            mu_.boundaryField()[i][j] = materialMixturePtr_->mu(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return mu_;
}


const volScalarField MNP::pv() const
{
    volScalarField pv_(IOobject(("pv_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("pv_"+particleName_), dimPressure, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        pv_.internalField()[i] = materialMixturePtr_->pv(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            pv_.boundaryField()[i][j] = materialMixturePtr_->pv(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return pv_;
}


const volScalarField MNP::hm() const
{
    volScalarField hm_(IOobject(("hm_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("hm_"+particleName_), dimEnthalpy, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        hm_.internalField()[i] = materialMixturePtr_->hm(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            hm_.boundaryField()[i][j] = materialMixturePtr_->hm(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return hm_;
}


const volScalarField MNP::Tm() const
{
    volScalarField Tm_(IOobject(("Tm_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tm_"+particleName_), dimTemperature, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tm_.internalField()[i] = materialMixturePtr_->Tm(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tm_.boundaryField()[i][j] = materialMixturePtr_->Tm(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return Tm_;
}


const volScalarField MNP::Tmr() const
{
    volScalarField Tmr_(IOobject(("Tmr_"+particleName_), T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar(("Tmr_"+particleName_), dimTemperature, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tmr_.internalField()[i] = materialMixturePtr_->Tmr(pi_, Ti_, MNPMoleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tmr_.boundaryField()[i][j] = materialMixturePtr_->Tmr(pij_, Tij_, MNPMoleFraction_);
	}
    }

    return Tmr_;
}


const volScalarField MNP::molten() const
{
    volScalarField molten_(IOobject("molten", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("molten", dimless, SMALL));
    return molten_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace particleModels
} // End namespace Foam

// ************************************************************************* //
