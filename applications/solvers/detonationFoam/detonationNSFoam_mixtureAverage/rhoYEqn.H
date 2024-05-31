tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        rhoPhi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

vk *= 0;
diffsh *= 0;
DiffError *= 0;
surfaceScalarField sumOfOtherSpeices = fvc::interpolate(DiffError) & mesh.Sf();

forAll(Y,i)
{
	volScalarField& Yi = Y[i];
	volScalarField& Xi = X[i];
	volScalarField& Diffi = Diff[i];
	DiffError += (rho*Diffi*fvc::grad(Yi) + Yi*rho*Diffi*fvc::grad(MW)/MW);
	if(thermalDiffusion)
	{
		if(Y[i].name() == "H")
		{
		   DiffError += (-rho*Diffi*TDRatio_H*fvc::grad(T)*Yi/(T*(Xi+SMALL)));
		}
		else if(Y[i].name() == "H2")
		{
		   DiffError +=  (-rho*Diffi*TDRatio_H2*fvc::grad(T)*Yi/(T*(Xi+SMALL)));
		}	
	}
}

reaction->correct();
Qdot = reaction->Qdot();
volScalarField Yt(0.0*Y[0]);
surfaceScalarField  phiUc = fvc::interpolate(DiffError) & mesh.Sf();

forAll(Y, i)
{
	if(i != inertIndex)
	{
		volScalarField& Yi = Y[i];
		volScalarField& Xi = X[i];
		volScalarField& Diffi = Diff[i];
		volScalarField  ThermalDiffRatio = TDRatio_H*0;
		
		if(thermalDiffusion)
		{
			if(Y[i].name() == "H")
			{
			   ThermalDiffRatio = TDRatio_H;
			}
			else if(Y[i].name() == "H2")
			{
			   ThermalDiffRatio = TDRatio_H2;
			}		   
			else
			{
			   ThermalDiffRatio *= 0;
			}
		}			
		
				for(label icell = 0; icell < rho.size(); icell++)
		{
			Hsi[icell] = composition.Hs(i, p[icell], T[icell]);
			HEi[icell] = composition.Es(i, p[icell], T[icell]);
		}

		volScalarField::Boundary& HsiBf = Hsi.boundaryFieldRef();
		volScalarField::Boundary& HEiBf = HEi.boundaryFieldRef();
		
		forAll(HsiBf, patchi)
		{
			scalarField& HsiPatch = HsiBf[patchi];
			scalarField& HEiPatch = HEiBf[patchi];
			const scalarField& pp = p.boundaryField()[patchi];
			const scalarField& Tp = T.boundaryField()[patchi];

			forAll(HsiPatch, facei)
			{
				HsiPatch[facei] = composition.Hs(i, pp[facei], Tp[facei]);
				HEiPatch[facei] = composition.Es(i, pp[facei], Tp[facei]);
			}
		}
		surfaceScalarField speicesFlux = fvc::interpolate(rho*Diffi*fvc::grad(Yi)+Yi*rho*Diffi*fvc::grad(MW)/MW+(-rho*Diff[i]*ThermalDiffRatio*fvc::grad(T)*Yi/(T*(Xi+SMALL)))-Yi*DiffError)& mesh.Sf();
		vk += fvc::interpolate(Hsi)*speicesFlux;
		sumOfOtherSpeices += speicesFlux;
		diffsh += fvc::laplacian(mixAlpha*HEi, Yi);

		fvScalarMatrix YEqn
		(
			fvm::ddt(rho, Yi)
		  + mvConvection->fvmDiv(rhoPhi, Yi)
		  - fvm::laplacian(rho*Diffi, Yi) 
		  + mvConvection->fvmDiv(phiUc,Yi)     		  
		  - fvc::laplacian(rho*Diffi*Yi/MW,MW)
		  - fvc::laplacian(-rho*Diffi*ThermalDiffRatio*Yi/(T*(Xi+SMALL)),T)
		  ==
			reaction->R(Yi)
		);

	    YEqn.relax();

	    YEqn.solve(mesh.solver("Yi"));
	    
	    Yi.max(0.0);
	    Yi.correctBoundaryConditions();
	    Yt += Yi;
	}
}

{
	label i = inertIndex;
	//volScalarField& Yi = Y[i];
	for(label icell = 0; icell < rho.size(); icell++)
	{
		Hsi[icell] = composition.Hs(i, p[icell], T[icell]);
		HEi[icell] = composition.Es(i, p[icell], T[icell]);
	}

	volScalarField::Boundary& HsiBf = Hsi.boundaryFieldRef();
	volScalarField::Boundary& HEiBf = HEi.boundaryFieldRef();
	
	forAll(HsiBf, patchi)
	{
		scalarField& HsiPatch = HsiBf[patchi];
		scalarField& HEiPatch = HEiBf[patchi];
		const scalarField& pp = p.boundaryField()[patchi];
		const scalarField& Tp = T.boundaryField()[patchi];

		forAll(HsiPatch, facei)
		{
			HsiPatch[facei] = composition.Hs(i, pp[facei], Tp[facei]);
			HEiPatch[facei] = composition.Es(i, pp[facei], Tp[facei]);
		}
	}

	vk += fvc::interpolate(Hsi)*-sumOfOtherSpeices;
	diffsh += fvc::laplacian(mixAlpha*HEi, 1-Yt);	
}

Y[inertIndex] = scalar(1) - Yt;
Y[inertIndex].max(0.0);
