/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.1                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//------------------------------- nanoPCM4Foam project -------------------------------//
//Author
    //Ehsan Golab, SUT. All rights reserved.
    //Ehsan1996Golab@gmail.com

//--------------------------------------------------------------------------------------//

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])

m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])

//dimensions of internal square [x, y, z]->[W, H, L]
m4_define(intSqWB, 1.0)
m4_define(intSqWT, 2.0)
m4_define(intSqH, 8.0)

m4_define(intSqWN, 10)
m4_define(intSqHN, 80)


//dimensions of internal cylinder
m4_define(intCyRB, 2.0)
m4_define(intCyRT, 3.0)
m4_define(intCyH, 8.0)


//dimensions of external cylinder
m4_define(extCyRB, 3.0)
m4_define(extCyRT, 4.0)
m4_define(extCyH, 8.0)

m4_define(extCyRN, 10)
m4_define(extCyHN, 80)


//dimensions of external square
m4_define(extSqWB, 10.0)
m4_define(extSqLB, 12.0)
m4_define(extSqWT, 12.0)
m4_define(extSqLT, 14.0)
m4_define(extSqH, 8.0)

m4_define(extSqWN, 10)
m4_define(extSqLN, 20)
m4_define(extSqHN, 80)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 1.0;

vertices
(
  //-Internal square
    //-Bottom of square
    (    calc(intSqWB/2)    calc(-1*intSqH/2)    calc(-1*intSqWB/2)    ) vlabel(intSqB0)
    (    calc(intSqWB/2)    calc(-1*intSqH/2)    calc(intSqWB/2)    ) vlabel(intSqB1)
    (    calc(-1*intSqWB/2)    calc(-1*intSqH/2)    calc(intSqWB/2)    ) vlabel(intSqB2)
    (    calc(-1*intSqWB/2)    calc(-1*intSqH/2)    calc(-1*intSqWB/2)    ) vlabel(intSqB3)
    //-Top of square
    (    calc(intSqWT/2)    calc(intSqH/2)    calc(-1*intSqWT/2)    ) vlabel(intSqT0)
    (    calc(intSqWT/2)    calc(intSqH/2)    calc(intSqWT/2)    ) vlabel(intSqT1)
    (    calc(-1*intSqWT/2)    calc(intSqH/2)    calc(intSqWT/2)    ) vlabel(intSqT2)
    (    calc(-1*intSqWT/2)    calc(intSqH/2)    calc(-1*intSqWT/2)    ) vlabel(intSqT3)


  //-Internal cylinder
    //-Bottom of internal cylinder
    (    calc(intCyRB*cos(rad(45.0)))    calc(-1*intCyH/2)    calc(-1*intCyRB*cos(rad(45.0)))    ) vlabel(intCyB0)
    (    calc(intCyRB*cos(rad(45.0)))    calc(-1*intCyH/2)    calc(intCyRB*cos(rad(45.0)))    ) vlabel(intCyB1)
    (    calc(-1*intCyRB*cos(rad(45.0)))    calc(-1*intCyH/2)    calc(intCyRB*cos(rad(45.0)))    ) vlabel(intCyB2)
    (    calc(-1*intCyRB*cos(rad(45.0)))    calc(-1*intCyH/2)    calc(-1*intCyRB*cos(rad(45.0)))    ) vlabel(intCyB3)
    //-Top of internal cylinder
    (    calc(intCyRT*cos(rad(45.0)))    calc(intCyH/2)    calc(-1*intCyRT*cos(rad(45.0)))    ) vlabel(intCyT0)
    (    calc(intCyRT*cos(rad(45.0)))    calc(intCyH/2)    calc(intCyRT*cos(rad(45.0)))    ) vlabel(intCyT1)
    (    calc(-1*intCyRT*cos(rad(45.0)))    calc(intCyH/2)    calc(intCyRT*cos(rad(45.0)))    ) vlabel(intCyT2)
    (    calc(-1*intCyRT*cos(rad(45.0)))    calc(intCyH/2)    calc(-1*intCyRT*cos(rad(45.0)))    ) vlabel(intCyT3)


  //-External cylinder
    //-Bottom of external cylinder
    (    calc(extCyRB*cos(rad(45.0)))    calc(-1*extCyH/2)    calc(-1*extCyRB*cos(rad(45.0)))    ) vlabel(extCyB0)
    (    calc(extCyRB*cos(rad(45.0)))    calc(-1*extCyH/2)    calc(extCyRB*cos(rad(45.0)))    ) vlabel(extCyB1)
    (    calc(-1*extCyRB*cos(rad(45.0)))    calc(-1*extCyH/2)    calc(extCyRB*cos(rad(45.0)))    ) vlabel(extCyB2)
    (    calc(-1*extCyRB*cos(rad(45.0)))    calc(-1*extCyH/2)    calc(-1*extCyRB*cos(rad(45.0)))    ) vlabel(extCyB3)
    //-Top of external cylinder
    (    calc(extCyRT*cos(rad(45.0)))    calc(extCyH/2)    calc(-1*extCyRT*cos(rad(45.0)))    ) vlabel(extCyT0)
    (    calc(extCyRT*cos(rad(45.0)))    calc(extCyH/2)    calc(extCyRT*cos(rad(45.0)))    ) vlabel(extCyT1)
    (    calc(-1*extCyRT*cos(rad(45.0)))    calc(extCyH/2)    calc(extCyRT*cos(rad(45.0)))    ) vlabel(extCyT2)
    (    calc(-1*extCyRT*cos(rad(45.0)))    calc(extCyH/2)    calc(-1*extCyRT*cos(rad(45.0)))    ) vlabel(extCyT3)


  //-External square
    //-Bottom of square
    (    calc(extSqWB/2)    calc(-1*extSqH/2)    calc(-1*extSqLB/2)    ) vlabel(extSqB0)
    (    calc(extSqWB/2)    calc(-1*extSqH/2)    calc(extSqLB/2)    ) vlabel(extSqB1)
    (    calc(-1*extSqWB/2)    calc(-1*extSqH/2)    calc(extSqLB/2)    ) vlabel(extSqB2)
    (    calc(-1*extSqWB/2)    calc(-1*extSqH/2)    calc(-1*extSqLB/2)    ) vlabel(extSqB3)
    //-Top of square
    (    calc(extSqWT/2)    calc(extSqH/2)    calc(-1*extSqLT/2)    ) vlabel(extSqT0)
    (    calc(extSqWT/2)    calc(extSqH/2)    calc(extSqLT/2)    ) vlabel(extSqT1)
    (    calc(-1*extSqWT/2)    calc(extSqH/2)    calc(extSqLT/2)    ) vlabel(extSqT2)
    (    calc(-1*extSqWT/2)    calc(extSqH/2)    calc(-1*extSqLT/2)    ) vlabel(extSqT3)
);

blocks
(
  //-Internal square
    hex (intSqB1 intSqB2 intSqT2 intSqT1 intSqB0 intSqB3 intSqT3 intSqT0) PCM (intSqWN intSqHN intSqWN) simpleGrading (1 1 1)

  //-Internal cylinder
    hex (intSqB0 intSqB3 intSqT3 intSqT0 intCyB0 intCyB3 intCyT3 intCyT0) PCM (intSqWN intSqHN intSqWN) simpleGrading (1 1 1)
    hex (intSqB3 intSqB2 intSqT2 intSqT3 intCyB3 intCyB2 intCyT2 intCyT3) PCM (intSqWN intSqHN intSqWN) simpleGrading (1 1 1)
    hex (intSqB2 intSqB1 intSqT1 intSqT2 intCyB2 intCyB1 intCyT1 intCyT2) PCM (intSqWN intSqHN intSqWN) simpleGrading (1 1 1)
    hex (intSqB1 intSqB0 intSqT0 intSqT1 intCyB1 intCyB0 intCyT0 intCyT1) PCM (intSqWN intSqHN intSqWN) simpleGrading (1 1 1)

  //-External cylinder
    hex (intCyB0 intCyB3 intCyT3 intCyT0 extCyB0 extCyB3 extCyT3 extCyT0) metalCover (extCyRN extCyHN extCyRN) simpleGrading (1 1 1)
    hex (intCyB3 intCyB2 intCyT2 intCyT3 extCyB3 extCyB2 extCyT2 extCyT3) metalCover (extCyRN extCyHN extCyRN) simpleGrading (1 1 1)
    hex (intCyB2 intCyB1 intCyT1 intCyT2 extCyB2 extCyB1 extCyT1 extCyT2) metalCover (extCyRN extCyHN extCyRN) simpleGrading (1 1 1)
    hex (intCyB1 intCyB0 intCyT0 intCyT1 extCyB1 extCyB0 extCyT0 extCyT1) metalCover (extCyRN extCyHN extCyRN) simpleGrading (1 1 1)

  //-External square
    hex (extCyB0 extCyB3 extCyT3 extCyT0 extSqB0 extSqB3 extSqT3 extSqT0) brick (extSqWN extSqHN extSqLN) simpleGrading (1 1 1)
    hex (extCyB3 extCyB2 extCyT2 extCyT3 extSqB3 extSqB2 extSqT2 extSqT3) brick (extSqWN extSqHN extSqLN) simpleGrading (1 1 1)
    hex (extCyB2 extCyB1 extCyT1 extCyT2 extSqB2 extSqB1 extSqT1 extSqT2) brick (extSqWN extSqHN extSqLN) simpleGrading (1 1 1)
    hex (extCyB1 extCyB0 extCyT0 extCyT1 extSqB1 extSqB0 extSqT0 extSqT1) brick (extSqWN extSqHN extSqLN) simpleGrading (1 1 1)
);

edges
(
  //-Internal cylinder
    //-Bottom of internal cylinder
    arc intCyB0 intCyB1 (    intCyRB    calc(-1*intCyH/2)    0    )
    arc intCyB1 intCyB2 (    0    calc(-1*intCyH/2)    intCyRB    )
    arc intCyB2 intCyB3 (    calc(-1*intCyRB)    calc(-1*intCyH/2)    0    )
    arc intCyB3 intCyB0 (    0    calc(-1*intCyH/2)    calc(-1*intCyRB)    )
    //-Top of internal cylinder
    arc intCyT0 intCyT1 (    intCyRT    calc(intCyH/2)    0    )
    arc intCyT1 intCyT2 (    0    calc(intCyH/2)    intCyRT    )
    arc intCyT2 intCyT3 (    calc(-1*intCyRT)    calc(intCyH/2)    0    )
    arc intCyT3 intCyT0 (    0    calc(intCyH/2)    calc(-1*intCyRT)    )

  //-External cylinder
    arc extCyB0 extCyB1 (    extCyRB    calc(-1*extCyH/2)    0    )
    arc extCyB1 extCyB2 (    0    calc(-1*extCyH/2)    extCyRB    )
    arc extCyB2 extCyB3 (    calc(-1*extCyRB)    calc(-1*extCyH/2)    0    )
    arc extCyB3 extCyB0 (    0    calc(-1*extCyH/2)    calc(-1*extCyRB)    )
    //-Top of external cylinder
    arc extCyT0 extCyT1 (    extCyRT    calc(extCyH/2)    0    )
    arc extCyT1 extCyT2 (    0    calc(extCyH/2)    extCyRT    )
    arc extCyT2 extCyT3 (    calc(-1*extCyRT)    calc(extCyH/2)    0    )
    arc extCyT3 extCyT0 (    0    calc(extCyH/2)    calc(-1*extCyRT)    )
);

boundary
(
    bottom
    {
        type wall;
        faces
        (
          //-Internal square
            (intSqB0 intSqB1 intSqB2 intSqB3)
          //-Internal cylinder
            (intSqB0 intCyB0 intCyB1 intSqB1)
            (intSqB3 intCyB3 intCyB0 intSqB0)
            (intSqB2 intCyB2 intCyB3 intSqB3)
            (intSqB1 intCyB1 intCyB2 intSqB2)
          //-External cylinder
            (intCyB0 extCyB0 extCyB1 intCyB1)
            (intCyB3 extCyB3 extCyB0 intCyB0)
            (intCyB2 extCyB2 extCyB3 intCyB3)
            (intCyB1 extCyB1 extCyB2 intCyB2)
          //-External square
            (extCyB0 extSqB0 extSqB1 extCyB1)
            (extCyB3 extSqB3 extSqB0 extCyB0)
            (extCyB2 extSqB2 extSqB3 extCyB3)
            (extCyB1 extSqB1 extSqB2 extCyB2)
        );
    }

    top
    {
        type wall;
        faces
        (
          //-Internal square
            (intSqT0 intSqT1 intSqT2 intSqT3)
          //-Internal cylinder
            (intSqT0 intCyT0 intCyT1 intSqT1)
            (intSqT3 intCyT3 intCyT0 intSqT0)
            (intSqT2 intCyT2 intCyT3 intSqT3)
            (intSqT1 intCyT1 intCyT2 intSqT2)
          //-External cylinder
            (intCyT0 extCyT0 extCyT1 intCyT1)
            (intCyT3 extCyT3 extCyT0 intCyT0)
            (intCyT2 extCyT2 extCyT3 intCyT3)
            (intCyT1 extCyT1 extCyT2 intCyT2)
          //-External square
            (extCyT0 extSqT0 extSqT1 extCyT1)
            (extCyT3 extSqT3 extSqT0 extCyT0)
            (extCyT2 extSqT2 extSqT3 extCyT3)
            (extCyT1 extSqT1 extSqT2 extCyT2)
        );
    }

    left
    {
        type wall;
        faces
        (
            (extSqB3 extSqB2 extSqT2 extSqT3)
        );
    }

    right
    {
        type wall;
        faces
        (
            (extSqB1 extSqB0 extSqT0 extSqT1)
        );
    }

    back
    {
        type wall;
        faces
        (
            (extSqB0 extSqB3 extSqT3 extSqT0)
        );
    }

    front
    {
        type wall;
        faces
        (
            (extSqB2 extSqB1 extSqT1 extSqT2)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
