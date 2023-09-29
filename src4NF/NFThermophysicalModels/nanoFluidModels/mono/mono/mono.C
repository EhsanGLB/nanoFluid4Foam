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
#include "mono.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nanoFluidModels
{
    defineTypeNameAndDebug(mono, 0);
    addToRunTimeSelectionTable(nanoFluidModel, mono, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mono::mono
(
    const dictionary& nanoFluidPropertiesDict,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& alphasPtr
)
:
    nanoFluidModel(nanoFluidPropertiesDict, U, p, T, alphasPtr),
    particleModelPtr_(particleModel::New(particlesProperties_[0].particleList(), nanoFluidPropertiesDict, U, p, T)),
    monoThermalConductivityModelPtr_(monoThermalConductivityModel::New(nanoFluidPropertiesDict, U, p, T, alphasPtr[0], baseFluidPtr_, particleModelPtr_)),
    monoDynamicViscosityModelPtr_(monoDynamicViscosityModel::New(nanoFluidPropertiesDict, U, p, T, alphasPtr[0], baseFluidPtr_, particleModelPtr_))
{
    if( particlesComponents_.size() != 1 )
    {
        FatalErrorIn("Stopped") << boldRed << "The nanofluid model is mono, so you should insert one component for particle." << reset << "\n" << exit(FatalError);
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const volScalarField mono::rho() const
{
    volScalarField rho_(IOobject("rhoNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhoNF", dimDensity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField rhonp_ = particlesProperties_[0].rho();

    rho_ = alpha() * rhonp_ + ( 1 - alpha() ) * rhobf_;

    return rho_;
}


const volScalarField mono::kappa() const
{
    volScalarField kappa_(IOobject("kappaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("kappaNF", dimThermalConductivity, SMALL));
    volScalarField kappabf_ = baseFluidPtr_->kappa();
    volScalarField kappanp_ = particlesProperties_[0].kappa();//change

    kappa_ = monoThermalConductivityModelPtr_->kappa();

    return kappa_;
}


const volScalarField mono::Cp() const
{
    volScalarField Cp_(IOobject("CpNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("CpNF", dimSpecificHeatCapacity, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField rhonp_ = particlesProperties_[0].rho();
    volScalarField Cpbf_ = baseFluidPtr_->Cp();
    volScalarField Cpnp_ = particlesProperties_[0].Cp();

    Cp_ = ( alpha() * rhonp_ * Cpnp_ + ( 1 - alpha() ) * rhobf_ * Cpbf_ ) / rho();

    return Cp_;
}


const volScalarField mono::beta() const
{
    volScalarField beta_(IOobject("betaNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("betaNF", dimThermalExpansion, SMALL));
    volScalarField rhobf_ = baseFluidPtr_->rho();
    volScalarField rhonp_ = particlesProperties_[0].rho();
    volScalarField betabf_ = baseFluidPtr_->beta();
    volScalarField betanp_ = particlesProperties_[0].beta();

    beta_ = (alpha() * rhonp_ * betanp_ + ( 1 - alpha() ) * rhobf_ * betabf_ ) / rho() ;

    return beta_;
}


const volScalarField mono::rhoR() const
{
    volScalarField rhoR_(IOobject("rhoRNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("rhoRNF", dimElectricalResistivity, SMALL));
    volScalarField rhoRbf_ = baseFluidPtr_->rhoR();
    volScalarField rhoRnp_ = particlesProperties_[0].rhoR();

    rhoR_ = rhoRbf_ / ( 1 + 3 * alpha() * ( rhoRbf_ / rhoRnp_ - 1 ) / ( ( rhoRbf_ / rhoRnp_ + 2 ) - ( rhoRbf_ / rhoRnp_ - 1 ) * alpha() ) );

    return rhoR_;
}


const volScalarField mono::mu() const
{
    volScalarField mu_(IOobject("muNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("muNF", dimDynamicViscosity, SMALL));

    mu_ = monoDynamicViscosityModelPtr_->mu();

    return mu_;
}


const volScalarField mono::pv() const
{
    volScalarField pv_(IOobject("pvNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("pvNF", dimPressure, SMALL));
    volScalarField pvbf_ = baseFluidPtr_->pv();
    volScalarField pvnp_ = particlesProperties_[0].pv();

    pv_ = pvbf_;

    return pv_;
}


const volScalarField mono::hm() const
{
    volScalarField hm_(IOobject("hmNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("hmNF", dimEnthalpy, SMALL));
    volScalarField hmbf_ = baseFluidPtr_->hm();
    volScalarField hmnp_ = particlesProperties_[0].hm();

    hm_ = ( 1 - alpha() ) * hmbf_;

    return hm_;
}


const volScalarField mono::Tm() const
{
    volScalarField Tm_(IOobject("TmNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("TmNF", dimTemperature, SMALL));
    volScalarField Tmbf_ = baseFluidPtr_->Tm();
    volScalarField Tmnp_ = particlesProperties_[0].Tm();

    Tm_ = Tmbf_;

    return Tm_;
}


const volScalarField mono::Tmr() const
{
    volScalarField Tmr_(IOobject("TmrNF", T_.time().timeName(), T_.db(), IOobject::NO_READ, IOobject::NO_WRITE), T_.mesh(), dimensionedScalar("TmrNF", dimTemperature, SMALL));
    volScalarField Tmrbf_ = baseFluidPtr_->Tmr();
    volScalarField Tmrnp_ = particlesProperties_[0].Tmr();

    Tmr_ = Tmrbf_;

    return Tmr_;
}


const volScalarField mono::alpha() const
{
    return alphasPtr_[0];
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void mono::print() const
{
    volScalarField Tmnp_ = particlesProperties_[0].Tm();
    volScalarField alpha_ = alpha();
    volScalarField Tm_ = Tm();
    volScalarField pv_ = pv();

    bool n1 = false;
    bool n2 = false;
    bool n3 = false;
    forAll(T_, i)
    {
        if ( T_[i] <= Tm_[i] )
        {
            n1 = true;
        }

        if ( p_[i] <= pv_[i] )
        {
            n2 = true;
        }

        if ( (alpha_[i] > 0.0) & (T_[i] >= Tmnp_[i]) )
        {
            n3 = true;
        }
    }

    if ( n1 == true )
    {
        cout << boldRed << "Warning: " << reset << "Flow temperature is lower than the Fluid's melting temperature." << "\n";
    }

    if ( n2 == true )
    {
        cout << boldRed << "Warning: " << reset << "Flow pressure is lower than the Fluid's vapour pressure (i.e. saturation state)." << "\n";
    }

    if ( n3 == true )
    {
        cout << boldRed << "Warning: " << reset << "Flow temperature is higher than the particle's melting temperature (i.e. particles are melting)." << "\n";
    }

    Info << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nanoFluidModels
} // End namespace Foam

// ************************************************************************* //
