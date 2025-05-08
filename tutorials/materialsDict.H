//------------------------------- nanoFluid4Foam project -------------------------------//
//Author
    //Ehsan Golab, SUT. All rights reserved.
    //Ehsan1996Golab@gmail.com

//--------------------------------------------------------------------------------------//

//- Templates
//- A material with const properties
//{
//Name    newMaterialConst    coeffs    molecularWeight[kg/kmol]    criticalTemperature[K]    criticalPressure[Pa]    criticalVolume[m^3/mol]    meltingPoint[K]    meltingRange[K]    meltingEnthalpy[J/kg]    density[kg/m^3]    ThermalConductivity[W/(m.K)]    specificHeatCapacity[J/(kg.K)]    thermalExpansion[1/K]    electricalResistivity[ohm.m]    dynamicViscosity[Pa.s]    vapourPressure[Pa];
//}

//- A material with non-const properties
//NameDict
//{
//    newMaterialFuncT    newMaterialFuncT;
//    constCoeffs    molecularWeight[kg/kmol]    criticalTemperature[K]    criticalPressure[Pa]    criticalVolume[m^3/mol]    meltingPoint[K]    meltingRange[K]    meltingEnthalpy[J/kg]
//    rhoCoeffs      a0    a1    a2    a3    a4    a5    a6    b    c;
//    KCoeffs        a0    a1    a2    a3    a4    a5    a6    b    c;
//    cpCoeffs       a0    a1    a2    a3    a4    a5    a6    b    c;
//    betaCoeffs     a0    a1    a2    a3    a4    a5    a6    b    c;
//    rhoRCoeffs     a0    a1    a2    a3    a4    a5    a6    b    c;
//    muCoeffs       a0    a1    a2    a3    a4    a5    a6    b    c;
//    pvCoeffs       a0    a1    a2    a3    a4    a5    a6    b    c;
//}
//- f(T) = ( a0 + a1*T^1 + a2*T^2 + a3*T^3 + a4*T^4 + a5*T^5 + a6*T^6 ) * exp( b * T ) + c;


//**************************************************************************# Metallic material #**************************************************************************//
Cu    newMaterialConst    coeffs    63.54    1.0e-23    1.0e-23    1.0e-23    1356.0    2.0    205382.0    8900.0    398.0    386.0    17.6e-6    17.1e-9    4.25e-3    1.0e5;
Ag    newMaterialConst    coeffs    107.86    1.0e-23    1.0e-23    1.0e-23    1235.0    2.0    3315408.0    10500.0    429.0    235.0    18.9e-6    6.3e6    1.0e-23    1.0e5;

//**************************************************************************# Ceramic material #**************************************************************************//
AL2O3    newMaterialConst    coeffs    101.96    1.0e-23    1.0e-23    1.0e-23    2327.0    1.0    1160000.0    3970.0    25.0    703.0    7.7e-6    5.0e23    45.0e-3    1.0e5;
TiO2    newMaterialConst    coeffs    79.8    1.0e-23    1.0e-23    1.0e-23    2116.0    2.0    11200.0    4230.0    8.4    692.0    9.0e-6    1.1e7    1e-10    1.0e5;
CuO    newMaterialConst    coeffs    79.54    1.0e-23    1.0e-23    1.0e-23    1599.0    2.0    124000.0    6315.0    76.0    540.0    8.5e-6    10.0    1.0e-23    1.0e5;
SiO2    newMaterialConst    coeffs    79.54    1.0e-23    1.0e-23    1.0e-23    1599.0    2.0    124000.0    6315.0    76.0    540.0    8.5e-6    10.0    1.0e-23    1.0e5;

//**************************************************************************# Organic and Carbon-based material #***********************************************************//
n_octadecane    newMaterialConst    coeffs    254.5    747.0   13.0e5    0.0039    301.35    1.5    243500.0    841.0    0.1505    2127.0    9.1e-4    1.0e10    3.86e-3    1.0e-23;
n_Hexadecane    newMaterialConst    coeffs    226.4    722.0   14.1e5    1.0e-23    292.0    1.0    237000.0    770.0    0.2    2215.0    7.0e-4    1.0e-23    3.03e-3    1.0e-23;
n_Docosane    newMaterialConst    coeffs    310.6    786.0   9.8e5    1.0e-23    317.35    1.0    251612.0    778.0    0.2    2380.0    8.0e-4    1.0e10    1.34e-3    1.0e-23;
lauricAcid    newMaterialConst    coeffs    200.3    743.0   18.7e5    1.0e-23    319.0    4.7    187210.0    912.5    0.15    2285.0    9.0e-5    1.0e10    6.11e-3    1.0e-23;
aceticAcid    newMaterialConst    coeffs    60.05    593.0   57.8e5    1.0e-23    290.0    3.0    184000.0    1160.0    0.158    2050.0    1.3e-3    0.2    1.22e-3    1.0e-23;
MWCNT    newMaterialConst    coeffs    12.01    1.0e-23    1.0e-23    1.0e-23    3823.0    2.0    1e-10    1600.0    3000.0    796.0    8.5e-6    1.0e-7    4.25e-3    1.0e5;


//**************************************************************************# Polymeric material #**************************************************************************//
polyurethane    newMaterialConst    coeffs    88.1    1.0e-23    1.0e-23    1.0e-23    448.0    2.0    1e-10    1125.0    0.4    1250.0    4.0e-5    2.5e-12    6.0e-3    1.0e-23;
PMMA    newMaterialConst    coeffs    270.0    1.0e-23    1.0e-23    1.0e-23    433.0    2.0    1e-10    1185.0    0.2    1466.0    4.4e-4    1.0e-19    6.0e-3    1.0e-23;

//**************************************************************************# Inorganic material #**************************************************************************//
air    newMaterialConst    coeffs    28.96    132.63    3.78e+6    1.0e-23    1.0e-23    1.0e-23    1.0e-23    1.161    0.0263    1007.1    3.34e-3    1.0e+10    1.86e-5    1.0e-23;
water    newMaterialConst    coeffs    18.015    647.13    2.2055e+7    0.05595    273.15    1.0    333606.0    996.5    0.6105    4181.0    2.6e-4    2.0e+7    8.5e-4    1.0e5;
blood    newMaterialConst    coeffs    1e-23    1e-23    1e-23    1e-23    1e-23    1e-23    1e-23    1060    0.45    3770.0    1e-23    1e-23    3.2e-3    1e-23;
EG    newMaterialConst    coeffs    62.1    720.0    8.2e6    1.0e-23    260.2    1.0    180354.0    1114.5    0.25    2386.0    2.0e-5    15.0e-4    1.6e-2    1.0e5;

waterFuncTDict
{
//273 < T < 647 K in liquid phase
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    18.015    647.13    2.2055e+7    0.05595    273.15    1.0    333606.0;
    rhoCoeffs      -1.75978216e+04    2.69778230e+02    -1.60962040e+00    5.06575412e-03    -8.87787757e-06    8.20798917e-09    -3.13009110e-12    0.0    0.0;
    KCoeffs        -1.37099193e+01    1.90089563e-01    -1.06896723e-03    3.24556679e-06    -5.56845688e-09    5.09400587e-12    -1.94028829e-15    0.0    0.0;
    cpCoeffs       1.67900503e+05    -2.49438291e+03    1.57014696e+01    -5.22747798e-02    9.70949933e-05    -9.54338265e-08    3.88350812e-11    0.0    0.0;
    betaCoeffs     1.68932021e+00    -2.53342939e-02    1.55968519e-04    -5.04797875e-07    9.06542164e-10    -8.57058466e-13    3.33530560e-16    0.0    0.0;
    rhoRCoeffs     2.0e+7    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       2.71770436e-01    -3.45583231e-03    1.82339752e-05    -5.09636846e-08    7.94842476e-11    -6.55464650e-14    2.23221324e-17    0.0    0.0;
    pvCoeffs       6.28324861e+04    -9.07569992e+02    5.32804589e+00    -1.62418304e-02    2.71640351e-05    -2.39875604e-08    9.19493560e-12    0.0    0.0;
}

iceFuncTDict
{
//147 < T < 273 K in solid phase
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    18.015    647.13    2.2055e+7    0.05595    273.15    1.0    333606.0;
    rhoCoeffs      9.34576949e+02    2.18254111e-02    -3.23992217e-04    -1.57861617e-08    1.15264685e-10    -3.20466045e-16    -2.94265166e-17    0.0    0.0;
    KCoeffs        3.06044881e+00    -3.55293343e-03    -6.24777093e-11    4.06680056e-13    -1.47626083e-15    2.83404765e-18    -2.24835856e-21    0.0    0.0;
    cpCoeffs       1.06008564e+04    -3.22145980e+02    4.32385428e+00    -3.02302076e-02    1.18523456e-04    -2.47000917e-07    2.13860840e-10    0.0    0.0;
    betaCoeffs     -2.33591197e-05    6.94064160e-07    2.41478477e-11    -2.38996952e-13    -2.14808275e-17    1.05645018e-19    -3.15371984e-23    0.0    0.0;
    rhoRCoeffs     2.0e+5    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    pvCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
}



therminolXP    newMaterialConst    coeffs    350.0    815.15    15.2e+5    0.00357    244.0    10.0    1e-23    874.0    0.1237    1850.0    8.92e-4    5.0e+12    39.2e-3    1.0e5;
therminolXPFuncTDict
{
//244 < T < 600 K in liquid phase
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    350.0    815.15    15.2e+5    0.00357    244.0    10.0    1e-23;
    rhoCoeffs      1.40039809e+03    -5.95640725e+00    3.40837846e-02    -1.13374384e-04    2.07200261e-07    -1.97628154e-10    7.66388048e-14    0.0    0.0;
    KCoeffs        1.14660354e-01    1.48243441e-04    -5.96129360e-07    7.53116691e-10    4.37444729e-14    -1.33258948e-15    1.00306582e-18    0.0    0.0;
    cpCoeffs       3.65020564e+03    -4.69274640e+01    3.30116054e-01    -1.08369126e-03    1.94402740e-06    -1.82333296e-09    6.99087699e-13    0.0    0.0;
    betaCoeffs     8.94e-4    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    rhoRCoeffs     5.0e12    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       362078298473237.2    0.0    0.0    0.0    0.0    0.0    0.0    -0.1300092602620023    0.00966732024088783;
    pvCoeffs       -1.00400015e+05    1.28897831e+03    -5.69661827e+00    7.04536192e-03    1.83409162e-05    -6.07803632e-08    4.84233983e-11    0.0    0.0;
}



therminolVP1    newMaterialConst    coeffs    166.0    722.0    33.1e5    0.003058104    285.0    5.0    4700    1056.0    0.1353    1575.0    1.0e-23    1.0e-23    3.28e-3    4.0;//300 K
therminolVP1FuncTDict
{
//285 < T < 693 K in liquid phase
    newMaterialFuncT    newMaterialFuncT;
    constCoeffs    166.0    722.0    33.1e5    0.003058104    285.0    5.0    4700;
    rhoCoeffs      1.4386e+03    -1.8711e+00    2.737e-03    -2.3793e-06    0.0    0.0    0.0    0.0    0.0;
    KCoeffs        0.14644    2.0353e-05    -1.9367e-07    1.0614e-11    0.0    0.0    0.0    0.0    0.0;
    cpCoeffs       2.125e+03    -11.017    0.049862    -7.7663e-05    4.394e-08    0.0    0.0    0.0    0.0;
    betaCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    rhoRCoeffs     1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
    muCoeffs       4.15    0.0    0.0    0.0    0.0    0.0    0.0    -2.363e-02    2.5847e-04;
    pvCoeffs       1.0e-23    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0;
}
