//*****************S. Zhong 2016.12.26******************//
  Info << "Read binaryDiff" << endl;
     IOdictionary binaryDiff
    (
        IOobject
        (
            "binaryDiff",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

const dimensionedScalar dUnit
(
   "dUnit",
   dimensionSet(0,0,0,0,0,0,0),
    1.0
);

boolList flags(0.5*Y.size()*(Y.size()+1), true);
forAll(groupSpecies,i)
{
	dictionary dDict;
	if(groupSpecies[i][0] != -1)
	{
		forAll(groupSpecies[i],h)
		{
			for(label j = i; j < groupSpecies.size();j ++)
			{
				if(groupSpecies[j][0] != -1)
				{
					forAll(groupSpecies[j],z)
					{
						bool foundDict = false;
						if(binaryDiff.found(Y[groupSpecies[i][h]].name()+"-"+Y[groupSpecies[j][z]].name()))
						{
							dDict = binaryDiff.subDict(Y[groupSpecies[i][h]].name()+"-"+Y[groupSpecies[j][z]].name());
							foundDict = true;
						}
						else if(binaryDiff.found(Y[groupSpecies[j][z]].name()+"-"+Y[groupSpecies[i][h]].name()))
						{
							dDict = binaryDiff.subDict(Y[groupSpecies[j][z]].name()+"-"+Y[groupSpecies[i][h]].name());
							foundDict = true;
						}
						else 
						{
							label xx = z*h;
							if(xx == (groupSpecies[j].size()*groupSpecies[i].size()))
							{
								Info << "binary diffsusion coefficient for " << endl;
								Info << "(" << endl;
								forAll(groupSpecies[i],h)
								{
									Info << Y[groupSpecies[i][h]].name() << endl;
								}
								Info << ")" << endl;
								Info << "(" << endl;
								forAll(groupSpecies[j],z)
								{
									Info << Y[groupSpecies[j][z]].name() << endl;
								}
								Info << ")" << endl;
								Info << " are not provided in constant/binaryDiff" <<	FatalErrorIn( "something wrong in your constant/binaryDiff") << abort(FatalError);	
							}						
						}
						
						if(foundDict)
						{
/* 							if(i != groupSpecies[j][0])
							{
								forAll(groupSpecies[j],z)
								{
										label m = max(i,groupSpecies[j][z]);
										label n = min(i,groupSpecies[j][z]);
										
										label k = Y.size()*n+m-0.5*n*(n+1);
										if(flags[k])
										{
											Diff1[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
											Diff2[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
											Diff3[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
											Diff4[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

											Diff1[k] = readScalar(dDict.lookup("Diff1")) * dUnit;
											Diff2[k] = readScalar(dDict.lookup("Diff2")) * dUnit;
											Diff3[k] = readScalar(dDict.lookup("Diff3")) * dUnit;
											Diff4[k] = readScalar(dDict.lookup("Diff4")) * dUnit;
											flags[k] = false;
										}	
								}
							} */
								forAll(groupSpecies[i],h)
								{
									forAll(groupSpecies[j],z)
									{
											label m = max(groupSpecies[i][h],groupSpecies[j][z]);
											label n = min(groupSpecies[i][h],groupSpecies[j][z]);
											
											label k = Y.size()*n+m-0.5*n*(n+1);
											if(flags[k])
											{
												Diff1[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
												Diff2[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
												Diff3[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
												Diff4[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

												Diff1[k] = readScalar(dDict.lookup("Diff1")) * dUnit;
												Diff2[k] = readScalar(dDict.lookup("Diff2")) * dUnit;
												Diff3[k] = readScalar(dDict.lookup("Diff3")) * dUnit;
												Diff4[k] = readScalar(dDict.lookup("Diff4")) * dUnit;
												flags[k] = false;
											}	
									}								
								}
							break;				
						}
					}
				}	
			}				
		}
	}
}


Info << "species binary diffusion log polynomial fitting coefficient" << endl;
Info << "(" << endl;
forAll(Y,i)
{
   for (label j = i; j < Y.size(); j++)
   {
		label m = max(i,j);
		label n = min(i,j);
		label k = Y.size()*n+m-0.5*n*(n+1);
		Info << Y[i].name() << "-" << Y[j].name() << tab << Diff1[k].value() << tab << Diff2[k].value() << tab << Diff3[k].value() << tab << Diff4[k].value() << endl;
	}
}
Info << ")" << endl;
/*
forAll(Y,i)
{
   for (label j = i; j < Y.size(); j++)
   {
		label m = max(i,j);
		label n = min(i,j);
		dictionary dDict;
		if(binaryDiff.found(Y[m].name()+"-"+Y[n].name()))
		{
			dDict = binaryDiff.subDict(Y[m].name()+"-"+Y[n].name());
		}
		else if(binaryDiff.found(Y[n].name()+"-"+Y[m].name()))
		{
			dDict = binaryDiff.subDict(Y[n].name()+"-"+Y[m].name());
		}
		else 
		{
			Info << "binary diffsusion coefficient for " << Y[m].name() << " and "  << Y[n].name() << " are not provided in constant/binaryDiff" << endl;
			FatalErrorIn("something in your constant" )
            << abort(FatalError);
		}
		//dictionary dDict(binaryDiff.subDict(Y[m].name()+"-"+Y[n].name()));
		
		// dictionary binaryDiffusionDict(YjDict.subDict(Y[j].name()));
		label k = Y.size()*n+m-0.5*n*(n+1);
		Diff1[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
		Diff2[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
		Diff3[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
		Diff4[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

		Diff1[k] = readScalar(dDict.lookup("Diff1")) * dUnit;
		Diff2[k] = readScalar(dDict.lookup("Diff2")) * dUnit;
		Diff3[k] = readScalar(dDict.lookup("Diff3")) * dUnit;
		Diff4[k] = readScalar(dDict.lookup("Diff4")) * dUnit;
   }
}
*/

