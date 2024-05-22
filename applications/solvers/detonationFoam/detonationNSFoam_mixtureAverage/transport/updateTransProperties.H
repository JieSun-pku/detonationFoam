mixMu *= 0;
mixLambdaA *= 0 ; // reset the number to zero;
mixLambdaB *= 0 ;

MW *= 0;
forAll(Y,i)
{
    MW += Y[i]/composition.W(i);
}
MW = 1/MW;
forAll(Y,i)
{
    X[i] = Y[i]*MW/composition.W(i);
}

volScalarField tmpT = log(T/TUnit);

if(differentialDiffusion)
{
	forAll(groupSpecies,i)
	{
		if(groupSpecies[i][0] != -1)
		{
			for(label j = i; j < groupSpecies.size();j ++)
			{
				if(groupSpecies[j][0] != -1)
				{
					label m = max(groupSpecies[i][0],groupSpecies[j][0]);
					label n = min(groupSpecies[i][0],groupSpecies[j][0]);	
					label k = Y.size()*n+m-0.5*n*(n+1);
					volScalarField tmpd = 1e-4 * exp(Diff1[k]+tmpT*(Diff2[k]+tmpT*(Diff3[k]+tmpT*Diff4[k]))) / (p / 101325 /pUnit);
					forAll(groupSpecies[i],h)
					{
						forAll(groupSpecies[j],z)
						{
								label m = max(groupSpecies[i][h],groupSpecies[j][z]);
								label n = min(groupSpecies[i][h],groupSpecies[j][z]);	
								label k = Y.size()*n+m-0.5*n*(n+1);
								d[k] = tmpd;
						}
					}
				}
			}
		}
	}
}

forAll(groupSpecies,i)
{
	if(groupSpecies[i][0] != -1)
	{
		Lambda[groupSpecies[i][0]] = lambdaUnit*1e-5*exp(Lambda1[groupSpecies[i][0]] + tmpT*(Lambda2[groupSpecies[i][0]]+ tmpT*(Lambda3[groupSpecies[i][0]]+tmpT*Lambda4[groupSpecies[i][0]])));
		Mu[groupSpecies[i][0]] = 0.1*MuUnit*exp(Mu1[groupSpecies[i][0]] + tmpT*(Mu2[groupSpecies[i][0]]+tmpT*(Mu3[groupSpecies[i][0]]+tmpT*Mu4[groupSpecies[i][0]])));		
		forAll(groupSpecies[i],j)
		{
			if(j > 0)
			{
				Lambda[groupSpecies[i][j]] = Lambda[groupSpecies[i][0]];
				Mu[groupSpecies[i][j]] = Mu[groupSpecies[i][0]];
			}
		}		
	}
}

forAll(Y,i)
{
	if(differentialDiffusion)
	{
		DiNumerator *= 0;
		DiDenominator *= 0;
		for (label j = 0; j< Y.size(); j++)
		{
			label m = max(i,j);
			label n = min(i,j);
			volScalarField& Xj = X[j];

			if(j != i)
			{
			 label k = Y.size()*n+m-0.5*n*(n+1);
			 DiNumerator += (Xj+1e-12) * composition.W(j);
			 DiDenominator += MW*(Xj+1e-12)/d[k];
			} 
		}
		Diff[i] = DiNumerator*DUnit/ DiDenominator; //*********** Di m2/s*******************//
	}
}
	
forAll(Y,i)
{
     volScalarField& Xi = X[i];
     volScalarField& Mui = Mu[i];
     yhi *= 0;
     forAll (Y,z)
	 {
		volScalarField& Xz = X[z];
		volScalarField& Muz = Mu[z];
		yhi += Xz*0.353553*pow(1.0+Foam::sqrt(Mui/Muz)*pow025(composition.W(z)/composition.W(i)),2)/Foam::sqrt(1.0+composition.W(i)/composition.W(z));
	 }
     mixMu += (Xi)*Mui/yhi;
	 mixLambdaA += Xi*Lambda[i];
     mixLambdaB += Xi/Lambda[i];
}

mixAlpha = 0.5*(mixLambdaA+(1/mixLambdaB))/thermo.Cp();

if(thermalDiffusion)
{
	TDRatio_H *= 0;
	TDRatio_H2 *= 0;
	volScalarField tmpTT = T/TUnit;
	forAll(Y,i)
	{
		if(Y[i].name() == "H")
		{
			forAll(Y,j)
			{
				if(j != i)
				{
					TDRatio_H += X[i]*X[j]*(TDiff1[j] + tmpTT*(TDiff2[j]+ tmpTT*(TDiff3[j]+tmpTT*TDiff4[j])));
				}
			}			
		}
		
		if(Y[i].name() == "H2")
		{
			forAll(Y,j)
			{
				if(j != i)
				{
					TDRatio_H += X[i]*X[j]*(TDiff1[j+Y.size()] + tmpTT*(TDiff2[j+Y.size()]+ tmpTT*(TDiff3[j+Y.size()]+tmpTT*TDiff4[j+Y.size()])));
				}
			}			
		}		
	}
}


if(differentialDiffusion)
{
	if(thermalDiffusion)
	{
		Info << "WITH differential diffusion and  thermalDiffusion"  << endl;
	}
	else
	{
		Info << "WITH differential diffusion"  << endl;
	}
}
else
{
    Info << "WITHOUT differential diffusion "  << endl;
}
Info << "update transport properties" << endl;

