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

#include "lumpedLIBBC4NF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedLIBBC4NF::lumpedLIBBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    state_("discharge"),
    rho_(0.0),
    Cp_(0.0),
    current_(0.0),
    resistance_(0.0),
    dEdT_(0.0),
    TName_("T"),
    kappaEffName_("kappaEff"),
    initTemp_(0.0),
    SOC_(0.0),
    capacity_(0.0),
    volume_(0.0),
    currentForm_("polynomial"),
    amplitude_(0),
    frequency_(0),
    real_(0),
    imaginary_(0)
{}


Foam::lumpedLIBBC4NF::lumpedLIBBC4NF
(
    const lumpedLIBBC4NF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    state_(ptf.state_),
    rho_(ptf.rho_),
    Cp_(ptf.Cp_),
    current_(ptf.current_),
    resistance_(ptf.resistance_),
    dEdT_(ptf.dEdT_),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    initTemp_(ptf.initTemp_),
    SOC_(ptf.SOC_),
    capacity_(ptf.capacity_),
    volume_(ptf.volume_),
    currentForm_(ptf.currentForm_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    real_(ptf.real_),
    imaginary_(ptf.imaginary_)
{}


Foam::lumpedLIBBC4NF::lumpedLIBBC4NF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    state_(dict.lookup("state")),
    rho_(readScalar(dict.lookup("rho"))),
    Cp_(readScalar(dict.lookup("Cp"))),
    current_(),
    resistance_(),
    dEdT_(readScalar(dict.lookup("dEdT"))),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaEffName_(dict.lookupOrDefault<word>("kappaEff", "kappaEff")),
    initTemp_(readScalar(dict.lookup("initTemp"))),
    SOC_(readScalar(dict.lookup("SOC"))),
    capacity_(readScalar(dict.lookup("capacity"))),
    volume_(readScalar(dict.lookup("volume"))),
    currentForm_(dict.lookup("currentForm")),
    amplitude_(),
    frequency_(),
    real_(),
    imaginary_()
{

    Istream& isresistance_ = dict.lookup("resistance");
    isresistance_.format(IOstream::ASCII);
    isresistance_ >> resistance_;

    if (currentForm_ == "polynomial")
    {
        Istream& iscurrent_ = dict.lookup("current");
        iscurrent_.format(IOstream::ASCII);
        iscurrent_ >> current_;
    }
    else
    {
        Istream& isamplitude_ = dict.lookup("amplitude");
        isamplitude_.format(IOstream::ASCII);
        isamplitude_ >> amplitude_;

        Istream& isfrequency_ = dict.lookup("frequency");
        isfrequency_.format(IOstream::ASCII);
        isfrequency_ >> frequency_;

        Istream& isreal_ = dict.lookup("real");
        isreal_.format(IOstream::ASCII);
        isreal_ >> real_;

        Istream& isimaginary_ = dict.lookup("imaginary");
        isimaginary_.format(IOstream::ASCII);
        isimaginary_ >> imaginary_;
    }

}


Foam::lumpedLIBBC4NF::lumpedLIBBC4NF
(
    const lumpedLIBBC4NF& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    state_(fcvpvf.state_),
    rho_(fcvpvf.rho_),
    Cp_(fcvpvf.Cp_),
    current_(fcvpvf.current_),
    resistance_(fcvpvf.resistance_),
    dEdT_(fcvpvf.dEdT_),
    TName_(fcvpvf.TName_),
    kappaEffName_(fcvpvf.kappaEffName_),
    initTemp_(fcvpvf.initTemp_),
    SOC_(fcvpvf.SOC_),
    capacity_(fcvpvf.capacity_),
    volume_(fcvpvf.volume_),
    currentForm_(fcvpvf.currentForm_),
    amplitude_(fcvpvf.amplitude_),
    frequency_(fcvpvf.frequency_),
    real_(fcvpvf.real_),
    imaginary_(fcvpvf.imaginary_)
{}


Foam::lumpedLIBBC4NF::lumpedLIBBC4NF
(
    const lumpedLIBBC4NF& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    state_(ptf.state_),
    rho_(ptf.rho_),
    Cp_(ptf.Cp_),
    current_(ptf.current_),
    resistance_(ptf.resistance_),
    dEdT_(ptf.dEdT_),
    TName_(ptf.TName_),
    kappaEffName_(ptf.kappaEffName_),
    initTemp_(ptf.initTemp_),
    SOC_(ptf.SOC_),
    capacity_(ptf.capacity_),
    volume_(ptf.volume_),
    currentForm_(ptf.currentForm_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    real_(ptf.real_),
    imaginary_(ptf.imaginary_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::lumpedLIBBC4NF::resistanceFunction(scalar LIBTemp)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < resistance_.size() )
    {
        f_ += resistance_[n_]*pow(LIBTemp, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::lumpedLIBBC4NF::currentPolynomialFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < current_.size() )
    {
        f_ += current_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::complex Foam::lumpedLIBBC4NF::currentWavyFunction(scalar t)
{
    scalar omega_ = 2.0*mathematicalConstant::pi*frequency_;
    int n_ = 1;
    complex It_(real_[0], imaginary_[0]);

    while( n_ < real_.size() )
    {
        It_ += complex(real_[n_]*cos(n_*omega_*t), imaginary_[n_]*sin(n_*omega_*t));
        n_ += 1;
    }

    return amplitude_*It_;
}


Foam::scalar Foam::lumpedLIBBC4NF::currentFunction(scalar t)
{
    scalar f_ = 0.0;
    if (currentForm_ == "polynomial")
    {
        f_ = currentPolynomialFunction(t);
    }
    else
    {
        f_ = currentWavyFunction(t).Re();
    }
    return f_;
}


void Foam::lumpedLIBBC4NF::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    //- determine the sign of currernt
    scalar sign_;
    if(state_ == "discharge"){sign_ = -1.0;}
    else{sign_ = 1.0;}

    //- get the currernt time and time step
    scalar t_ = this->db().time().value();//this->db().time().timeOutputValue();
    scalar deltaT_ = this->db().time().deltaT().value();

    //- compute the SOC
    scalar SOCCurr_ = SOC_ + sign_*(100/(3600*capacity_))*currentFunction(t_)*deltaT_;
    if( SOCCurr_ <= 0.0 )
    {
        FatalErrorIn("finished") << exit(FatalError);
    }
    SOC_ = SOCCurr_;

    //- compute the LIB temperature by lumped model
    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappaEffp = lookupPatchField<volScalarField, scalar>(kappaEffName_);
    scalarField gradTp_ = Tp.snGrad();

    scalar QOhmic_ = resistanceFunction(initTemp_)*currentFunction(t_)*currentFunction(t_);
    scalar QEntropic_ = -1.0*sign_*currentFunction(t_)*initTemp_*dEdT_;
    scalar Qb_ = QOhmic_ + QEntropic_;

    scalar QLoss_ = 0.0;
    forAll(Tp, i)
    {
        QLoss_ += -kappaEffp[i]*gradTp_[i]*patch().magSf()[i];
    }

    if(t_<10){QLoss_ = 0.0;}
    scalar Tb_ = initTemp_ + deltaT_*( (Qb_+QLoss_) / (rho_*Cp_*volume_) );
    initTemp_ = Tb_;
 
    scalarField::operator=(Tb_);
    fixedValueFvPatchScalarField::updateCoeffs();

    Info << "LIB Temperature = " << Tb_ << endl;
    Info << "SOC = " << SOC_ << endl;
    Info << "QLoss = " << QLoss_ << endl;
    Info << "R = " << resistanceFunction(initTemp_) << endl;
}


// Write
void Foam::lumpedLIBBC4NF::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("state") << state_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cp") << Cp_ << token::END_STATEMENT << nl;
    os.writeKeyword("resistance") << resistance_ << token::END_STATEMENT << nl;
    os.writeKeyword("dEdT") << dEdT_ << token::END_STATEMENT << nl;
    os.writeKeyword("initTemp") << initTemp_ << token::END_STATEMENT << nl;
    os.writeKeyword("SOC") << SOC_ << token::END_STATEMENT << nl;
    os.writeKeyword("capacity") << capacity_ << token::END_STATEMENT << nl;
    os.writeKeyword("volume") << volume_ << token::END_STATEMENT << nl;

    os.writeKeyword("currentForm") << currentForm_ << token::END_STATEMENT << nl;
    if (currentForm_ == "polynomial")
    {
        os.writeKeyword("current") << current_ << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
        os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
        os.writeKeyword("real") << real_ << token::END_STATEMENT << nl;
        os.writeKeyword("imaginary") << imaginary_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField(fvPatchScalarField, lumpedLIBBC4NF);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

