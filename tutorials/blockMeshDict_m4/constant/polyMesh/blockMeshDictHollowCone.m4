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

//------------------------------- nanoFluid4Foam project -------------------------------//
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

//variables
m4_define(R1, 21)
m4_define(R2T, 36)
m4_define(R2B, 26)
m4_define(length, 100)
m4_define(numberOfCell_x, 20)
m4_define(numberOfCell_y, 20)
m4_define(numberOfCell_z, 50)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 0.001;

vertices
(
    //-Back
    (    calc(R1*cos(rad(0.0)))    calc(R1*sin(rad(0.0)))    calc(-1*length/2)    ) vlabel(B0)
    (    calc(R1*cos(rad(90.0)))    calc(R1*sin(rad(90.0)))    calc(-1*length/2)    ) vlabel(B1)
    (    calc(R1*cos(rad(180.0)))    calc(R1*sin(rad(180.0)))    calc(-1*length/2)    ) vlabel(B2)
    (    calc(R1*cos(rad(270.0)))    calc(R1*sin(rad(270.0)))    calc(-1*length/2) ) vlabel(B3)
    (    calc(R2B*cos(rad(0.0)))    calc(R2B*sin(rad(0.0)))    calc(-1*length/2)    ) vlabel(B4)
    (    calc(R2B*cos(rad(90.0)))    calc(R2B*sin(rad(90.0)))    calc(-1*length/2)    ) vlabel(B5)
    (    calc(R2B*cos(rad(180.0)))    calc(R2B*sin(rad(180.0)))    calc(-1*length/2)    ) vlabel(B6)
    (    calc(R2B*cos(rad(270.0)))    calc(R2B*sin(rad(270.0)))    calc(-1*length/2)    ) vlabel(B7)

    //-Front
    (    calc(R1*cos(rad(0.0)))    calc(R1*sin(rad(0.0)))    calc(length/2)    ) vlabel(F0)
    (    calc(R1*cos(rad(90.0)))    calc(R1*sin(rad(90.0)))    calc(length/2)    ) vlabel(F1)
    (    calc(R1*cos(rad(180.0)))    calc(R1*sin(rad(180.0)))    calc(length/2)    ) vlabel(F2)
    (    calc(R1*cos(rad(270.0)))    calc(R1*sin(rad(270.0)))    calc(length/2) ) vlabel(F3)
    (    calc(R2T*cos(rad(0.0)))    calc(R2T*sin(rad(0.0)))    calc(length/2)    ) vlabel(F4)
    (    calc(R2T*cos(rad(90.0)))    calc(R2T*sin(rad(90.0)))    calc(length/2)    ) vlabel(F5)
    (    calc(R2T*cos(rad(180.0)))    calc(R2T*sin(rad(180.0)))    calc(length/2)    ) vlabel(F6)
    (    calc(R2T*cos(rad(270.0)))    calc(R2T*sin(rad(270.0)))    calc(length/2)    ) vlabel(F7)
);

blocks
(
    hex (B0 B4 B5 B1 F0 F4 F5 F1) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
    hex (B1 B5 B6 B2 F1 F5 F6 F2) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
    hex (B2 B6 B7 B3 F2 F6 F7 F3) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
    hex (B3 B7 B4 B0 F3 F7 F4 F0) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
);

edges
(
    //-Back
    arc B0 B1 (    calc(R1*cos(rad(45.0)))    calc(R1*sin(rad(45.0)))    calc(-1*length/2)    )
    arc B1 B2 (    calc(R1*cos(rad(135.0)))    calc(R1*sin(rad(135.0)))    calc(-1*length/2)    )
    arc B2 B3 (    calc(R1*cos(rad(225.0)))    calc(R1*sin(rad(225.0)))    calc(-1*length/2)    )
    arc B3 B0 (    calc(R1*cos(rad(315.0)))    calc(R1*sin(rad(315.0)))    calc(-1*length/2)    )
    arc B4 B5 (    calc(R2B*cos(rad(45.0)))    calc(R2B*sin(rad(45.0)))    calc(-1*length/2)    )
    arc B5 B6 (    calc(R2B*cos(rad(135.0)))    calc(R2B*sin(rad(135.0)))    calc(-1*length/2)    )
    arc B6 B7 (    calc(R2B*cos(rad(225.0)))    calc(R2B*sin(rad(225.0)))    calc(-1*length/2)    )
    arc B7 B4 (    calc(R2B*cos(rad(315.0)))    calc(R2B*sin(rad(315.0)))    calc(-1*length/2)    )
    //-Front
    arc F0 F1 (    calc(R1*cos(rad(45.0)))    calc(R1*sin(rad(45.0)))    calc(length/2)    )
    arc F1 F2 (    calc(R1*cos(rad(135.0)))    calc(R1*sin(rad(135.0)))    calc(length/2)    )
    arc F2 F3 (    calc(R1*cos(rad(225.0)))    calc(R1*sin(rad(225.0)))    calc(length/2)    )
    arc F3 F0 (    calc(R1*cos(rad(315.0)))    calc(R1*sin(rad(315.0)))    calc(length/2)    )
    arc F4 F5 (    calc(R2T*cos(rad(45.0)))    calc(R2T*sin(rad(45.0)))    calc(length/2)    )
    arc F5 F6 (    calc(R2T*cos(rad(135.0)))    calc(R2T*sin(rad(135.0)))    calc(length/2)    )
    arc F6 F7 (    calc(R2T*cos(rad(225.0)))    calc(R2T*sin(rad(225.0)))    calc(length/2)    )
    arc F7 F4 (    calc(R2T*cos(rad(315.0)))    calc(R2T*sin(rad(315.0)))    calc(length/2)    )
);

boundary
(
    innerCylinder
    {
        type wall;
        faces
        (
            (B1 B0 F0 F1)
            (B2 B1 F1 F2)
            (B3 B2 F2 F3)
            (B0 B3 F3 F0)
        );
    }

    outerCylinder
    {
        type wall;
        faces
        (
            (B4 B5 F5 F4)
            (B5 B6 F6 F5)
            (B6 B7 F7 F6)
            (B7 B4 F4 F7)
        );
    }

    bottom
    {
        type wall;
        faces
        (
            (B0 B1 B5 B4)
            (B1 B2 B6 B5)
            (B2 B3 B7 B6)
            (B3 B0 B4 B7)
        );
    }

    top
    {
        type wall;
        faces
        (
            (F0 F1 F5 F4)
            (F1 F2 F6 F5)
            (F2 F3 F7 F6)
            (F3 F0 F4 F7)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
