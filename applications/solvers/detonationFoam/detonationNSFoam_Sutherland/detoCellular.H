volScalarField magGradrho("magGradrho",mag(fvc::grad(rho)));  

for (label cellI=0; cellI<mesh.C().size(); cellI++)
{
    if(maxp[cellI]<p[cellI])  maxp[cellI] = p[cellI] ;            
}

maxp.writeOpt() = IOobject::AUTO_WRITE;
magGradrho.writeOpt() = IOobject::AUTO_WRITE;
