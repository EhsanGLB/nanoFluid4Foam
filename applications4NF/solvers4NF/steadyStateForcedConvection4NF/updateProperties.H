rho = nanoFluidModelPtr->rho();
kappa = nanoFluidModelPtr->kappa();
Cp = nanoFluidModelPtr->Cp();
mu = nanoFluidModelPtr->mu();
nanoFluidModelPtr->print();

forAll(moltensPtr, i)
{
    moltensPtr[i] = particlesProperties[i].molten();
}

/*
beta = nanoFluidModelPtr->beta();
rhoR = nanoFluidModelPtr->rhoR();
pv = nanoFluidModelPtr->pv();
hm = nanoFluidModelPtr->hm();
Tm = nanoFluidModelPtr->Tm();
Tmr = nanoFluidModelPtr->Tmr();
*/
