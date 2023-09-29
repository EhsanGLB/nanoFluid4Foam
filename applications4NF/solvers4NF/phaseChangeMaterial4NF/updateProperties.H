rho = nanoFluidModelPtr->rho();
kappa = nanoFluidModelPtr->kappa();
Cp = nanoFluidModelPtr->Cp();
mu = nanoFluidModelPtr->mu();
beta = nanoFluidModelPtr->beta();
hm = nanoFluidModelPtr->hm();
Tm = nanoFluidModelPtr->Tm();
Tmr = nanoFluidModelPtr->Tmr();
nanoFluidModelPtr->print();

forAll(moltensPtr, i)
{
    moltensPtr[i] = particlesProperties[i].molten();
}

gamma = 0.5 * erf( 4.0 * ( T - Tm ) / Tmr ) + scalar(0.5);


