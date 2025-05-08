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

#include "materialMixture.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::materialMixture::TrMax = 0.999;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materialMixture::materialMixture
(
    const wordList& components,
    const dictionary& propertiesSubDict
)
:
    components_(components),
    properties_(components_.size())
{

    forAll(components_, i)
    {
        const dictionary* subDictPtr = propertiesSubDict.subDictPtr(components_[i] + "Dict");

        if (subDictPtr)
        {
            properties_.set
            (
                i,
                material::New(subDictPtr)
            );
        }

        else
        {
            properties_.set
            (
                i,
                material::New(propertiesSubDict.lookup(components_[i]))
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::materialMixture> Foam::materialMixture::New
(
    const wordList& components,
    const dictionary& propertiesSubDict
)
{
    return autoPtr<materialMixture>(new materialMixture(components, propertiesSubDict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::materialMixture::Tc
(
    const scalarField& x
) const
{

    scalar vTc = 0.0;
    scalar vc = 0.0;

    forAll(properties_, i)
    {
        scalar x1 = x[i]*properties_[i].Vc();
        vc += x1;
        vTc += x1*properties_[i].Tc();
    }

    return vTc/vc;
}


Foam::scalar Foam::materialMixture::W
(
    const scalarField& x
) const
{
    scalar W = 0.0;
    forAll(properties_, i)
    {
        W += x[i]*properties_[i].W();
    }

    return W;
}


Foam::scalarField Foam::materialMixture::Y
(
    const scalarField& X
) const
{
    scalarField Y(X/W(X));

    forAll(Y, i)
    {
        Y[i] *= properties_[i].W();
    }

    return Y;
}


Foam::scalarField Foam::materialMixture::X
(
    const scalarField& Y
) const
{
    scalarField X(Y.size());
    scalar Winv = 0.0;
    forAll(X, i)
    {
        Winv += Y[i]/properties_[i].W();
        X[i] = Y[i]/properties_[i].W();
    }

    return X/Winv;
}


Foam::scalar Foam::materialMixture::rho
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar v = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            scalar rho = SMALL + properties_[i].rho(Ti);
            v += x[i]*properties_[i].W()/rho;
        }
    }

    return W(x)/v;
}


Foam::scalar Foam::materialMixture::K
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // calculate superficial volume fractions phii
    scalarField phii(x.size(), 0.0);
    scalar pSum = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        scalar Vi = properties_[i].W()/properties_[i].rho(Ti);
        phii[i] = x[i]*Vi;
        pSum += phii[i];
    }

    forAll(phii, i)
    {
        phii[i] /= pSum;
    }

    scalar K = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        forAll(properties_, j)
        {
            scalar Tj = min(TrMax*properties_[j].Tc(), T);

            scalar Kij =
                2.0
               /(
                    1.0/properties_[i].K(Ti)
                  + 1.0/properties_[j].K(Tj)
                );
            K += phii[i]*phii[j]*Kij;
        }
    }

    return K;
}


Foam::scalar Foam::materialMixture::cp
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar cp = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            cp += x[i]*properties_[i].cp(Ti)*properties_[i].W();
        }
    }

    return cp/W(x);
}


Foam::scalar Foam::materialMixture::beta
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar beta = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            beta += x[i]*properties_[i].beta(Ti)*properties_[i].W();
        }
    }

    return beta/W(x);
}


Foam::scalar Foam::materialMixture::rhoR
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar rhoR = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            rhoR += x[i]*properties_[i].rhoR(Ti)*properties_[i].W();
        }
    }

    return rhoR/W(x);
}


Foam::scalar Foam::materialMixture::mu
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar mu = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            mu += x[i]*log(properties_[i].mu(Ti));
        }
    }

    return exp(mu);
}


Foam::scalar Foam::materialMixture::pv
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar pv = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            pv += x[i]*properties_[i].pv(Ti)*properties_[i].W();
        }
    }

    return pv/W(x);
}


Foam::scalar Foam::materialMixture::hm
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar hm = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            hm += x[i]*properties_[i].hm()*properties_[i].W();
        }
    }

    return hm/W(x);
}


Foam::scalar Foam::materialMixture::sm
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar sm = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            sm += x[i]*properties_[i].sm()*properties_[i].W();
        }
    }

    return sm/W(x);
}


Foam::scalar Foam::materialMixture::Tm
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    return hm(p, T, x)/sm(p, T, x);
}

Foam::scalar Foam::materialMixture::Tmr
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    return properties_[0].Tmr();
}

// ************************************************************************* //
