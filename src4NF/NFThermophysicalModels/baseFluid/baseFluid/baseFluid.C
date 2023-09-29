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

#include "volFields.H"

#include "NFDimensionSet.H"//-nanoFluid4Foam
#include "baseFluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<baseFluid> baseFluid::New
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
{
    return autoPtr<baseFluid>(new baseFluid(nanoFluidPropertiesDict, U, p, T));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

baseFluid::baseFluid
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T
)
:
    nanoFluidPropertiesDict_(nanoFluidPropertiesDict),
    U_(U),
    p_(p),
    T_(T),
    baseFluidComponents_(nanoFluidPropertiesDict.lookup("baseFluid")),
    moleFraction_(nanoFluidPropertiesDict.lookup("moleFraction")),
    materialMixturePtr_(materialMixture::New(baseFluidComponents_, nanoFluidPropertiesDict)),
    baseFluidDynamicViscosityModelPtr_(baseFluidDynamicViscosityModel::New(nanoFluidPropertiesDict, U, p, T))
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField baseFluid::rho() const
{
    volScalarField rho_(IOobject("rhoBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhoBF", dimDensity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rho_.internalField()[i] = materialMixturePtr_->rho(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rho_.boundaryField()[i][j] = materialMixturePtr_->rho(pij_, Tij_, moleFraction_);
	}
    }

    return rho_;
}


const volScalarField baseFluid::kappa() const
{
    volScalarField kappa_(IOobject("kappaBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaBF", dimThermalConductivity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        kappa_.internalField()[i] = materialMixturePtr_->K(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            kappa_.boundaryField()[i][j] = materialMixturePtr_->K(pij_, Tij_, moleFraction_);
	}
    }

    return kappa_;
}


const volScalarField baseFluid::Cp() const
{
    volScalarField Cp_(IOobject("CpBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("CpBF", dimSpecificHeatCapacity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Cp_.internalField()[i] = materialMixturePtr_->cp(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Cp_.boundaryField()[i][j] = materialMixturePtr_->cp(pij_, Tij_, moleFraction_);
	}
    }

    return Cp_;
}


const volScalarField baseFluid::beta() const
{
    volScalarField beta_(IOobject("betaBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("betaBF", dimThermalExpansion, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        beta_.internalField()[i] = materialMixturePtr_->beta(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            beta_.boundaryField()[i][j] = materialMixturePtr_->beta(pij_, Tij_, moleFraction_);
	}
    }

    return beta_;
}


const volScalarField baseFluid::rhoR() const
{
    volScalarField rhoR_(IOobject("rhoRBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhoRBF", dimElectricalResistivity, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        rhoR_.internalField()[i] = materialMixturePtr_->rhoR(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            rhoR_.boundaryField()[i][j] = materialMixturePtr_->rhoR(pij_, Tij_, moleFraction_);
	}
    }

    return rhoR_;
}


const volScalarField baseFluid::mu() const
{
    volScalarField mu_(IOobject("muBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muBF", dimDynamicViscosity, SMALL));
    mu_ = baseFluidDynamicViscosityModelPtr_->mu();

    return mu_;
}


const volScalarField baseFluid::pv() const
{
    volScalarField pv_(IOobject("pvBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("pvBF", dimPressure, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        pv_.internalField()[i] = materialMixturePtr_->pv(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            pv_.boundaryField()[i][j] = materialMixturePtr_->pv(pij_, Tij_, moleFraction_);
	}
    }

    return pv_;
}


const volScalarField baseFluid::hm() const
{
    volScalarField hm_(IOobject("hmBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("hmBF", dimEnthalpy, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        hm_.internalField()[i] = materialMixturePtr_->hm(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            hm_.boundaryField()[i][j] = materialMixturePtr_->hm(pij_, Tij_, moleFraction_);
	}
    }

    return hm_;
}


const volScalarField baseFluid::Tm() const
{
    volScalarField Tm_(IOobject("TmBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("TmBF", dimTemperature, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tm_.internalField()[i] = materialMixturePtr_->Tm(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tm_.boundaryField()[i][j] = materialMixturePtr_->Tm(pij_, Tij_, moleFraction_);
	}
    }

    return Tm_;
}


const volScalarField baseFluid::Tmr() const
{
    volScalarField Tmr_(IOobject("TmrBF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("TmrBF", dimTemperature, SMALL));
    forAll(T_.internalField(), i)
    {
        scalar pi_ = p_.internalField()[i];
        scalar Ti_ = T_.internalField()[i];
        Tmr_.internalField()[i] = materialMixturePtr_->Tmr(pi_, Ti_, moleFraction_);
    }

    forAll(T_.boundaryField(), i)
    {
	forAll(T_.boundaryField()[i], j)
	{
            scalar pij_ = p_.boundaryField()[i][j];
            scalar Tij_ = T_.boundaryField()[i][j];
            Tmr_.boundaryField()[i][j] = materialMixturePtr_->Tmr(pij_, Tij_, moleFraction_);
	}
    }

    return Tmr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
