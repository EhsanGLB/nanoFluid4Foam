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

//**************************************************************************# Ceramic material #**************************************************************************//
AL2O3    newMaterialConst    coeffs    101.96    1.0e-23    1.0e-23    1.0e-23    2327.0    1.0    1160000.0    3970.0    40.0    765.0    24e-6    5.0e23    45.0e-3    1.0e5;


//**************************************************************************# Organic and Carbon-based material #***********************************************************//


//**************************************************************************# Polymeric material #**************************************************************************//


//**************************************************************************# Inorganic material #**************************************************************************//
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

