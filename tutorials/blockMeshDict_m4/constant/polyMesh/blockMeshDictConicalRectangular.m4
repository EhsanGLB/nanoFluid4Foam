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

//dimensions [x, y, z]->[W, H, L]
//width and length of bottom
m4_define(WB, 1.0)
m4_define(LB, 2.0)
//width and length of top
m4_define(WT, 2.0)
m4_define(LT, 4.0)
//height
m4_define(H, 4.0)

// number of cell
m4_define(WN, 10)
m4_define(LN, 20)
m4_define(HN, 40)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 1.0;

vertices
(
    //-Bottom
    (    calc(WB/2)    calc(-1*H/2)    calc(-1*LB/2)    ) vlabel(B0)
    (    calc(WB/2)    calc(-1*H/2)    calc(LB/2)    ) vlabel(B1)
    (    calc(-1*WB/2)    calc(-1*H/2)    calc(LB/2)    ) vlabel(B2)
    (    calc(-1*WB/2)    calc(-1*H/2)    calc(-1*LB/2)    ) vlabel(B3)

    //-Top
    (    calc(WT/2)    calc(H/2)    calc(-1*LT/2)    ) vlabel(T0)
    (    calc(WT/2)    calc(H/2)    calc(LT/2)    ) vlabel(T1)
    (    calc(-1*WT/2)    calc(H/2)    calc(LT/2)    ) vlabel(T2)
    (    calc(-1*WT/2)    calc(H/2)    calc(-1*LT/2)    ) vlabel(T3)
);

blocks
(
    hex (B1 B2 T2 T1 B0 B3 T3 T0) (WN HN LN) simpleGrading (1 1 1)
);

edges
(

);

boundary
(
    bottom
    {
        type wall;
        faces
        (
            (B0 B1 B2 B3)
        );
    }

    top
    {
        type wall;
        faces
        (
            (T0 T3 T2 T1)
        );
    }

    left
    {
        type wall;
        faces
        (
            (B3 B2 T2 T3)
        );
    }

    right
    {
        type wall;
        faces
        (
            (B1 B0 T0 T1)
        );
    }

    back
    {
        type wall;
        faces
        (
            (B0 B3 T3 T0)
        );
    }

    front
    {
        type wall;
        faces
        (
            (B2 B1 T1 T2)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
