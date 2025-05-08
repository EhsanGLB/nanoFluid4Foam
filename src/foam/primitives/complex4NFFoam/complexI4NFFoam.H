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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline unsigned long int factorial(unsigned long int n)
{
    if (n == 0 || n == 1)
        return 1;
    return n * factorial(n - 1);
}


inline complex exp(const complex& c)
{
    return Foam::exp(c.Re())*complex(Foam::cos(c.Im()), Foam::sin(c.Im()));
}


inline complex log(const complex& c)
{
    complex log_;

    if (c.Re() == 0.0)
    {
        log_ = complex(0.5*Foam::log(magSqr(c)), Foam::atan(GREAT));
    }

    else
    {
        log_ = complex(0.5*Foam::log(magSqr(c)), Foam::atan(c.Im()/c.Re()));
    }

    return log_;
}


inline complex pow(const complex& c, const scalar s)
{
    complex pow_;

    if (c == complex(0, 0))
    {
        pow_ = complex(0.0, 0.0);
    }

    else
    {
        pow_ = exp(s * log(c));
    }

    return pow_;
}


inline complex sqrt(const complex& c)
{
    scalar re_ = Foam::sqrt((mag(c) + c.Re())/2.0);
    scalar im_ = sign(c.Im())*Foam::sqrt((mag(c) - c.Re())/2.0);

    return complex(re_, im_);
}


inline complex sin(const complex& c)
{
    complex i_(0.0, 1.0);

    return (exp(i_*c)-exp(-1.0*i_*c))/(2.0*i_);
}


inline complex cos(const complex& c)
{
    complex i_(0.0, 1.0);

    return (exp(i_*c)+exp(-1.0*i_*c))/(2.0);
}


inline complex J(const int& n, const complex& c)
{
    complex J_;

    if (c == complex(0, 0))
    {
        if(n == 0)
        {
            J_ = complex(1.0, 1.0);
        }
        else
        {
            J_ = complex(0.0, 0.0);
        }
    }

    else
    {
        int k_ = 0;

        while(k_ < 12)
        {
            J_ += (pow(-1.0, k_)/(factorial(k_)*factorial(k_+n)))*pow((c/2), 2.0*k_+n);
            k_ += 1;
        }
    }

    return J_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
