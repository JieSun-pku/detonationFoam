//********************S. Zhong 2017.3.14**************//
  Info << "Read speciesLambda" << endl;

    IOdictionary speciesLambda
    (
        IOobject
        (
            "speciesLambda",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

//******in order to give unit to the polynomial coefficient***//
const dimensionedScalar speciesLambdaUnit
(
   "speciesLambdaUnit",
   dimensionSet(0,0,0,0,0,0,0),
    1.0
);


forAll(groupSpecies,i)
{
	dictionary lambdaDict;
	if(groupSpecies[i][0] != -1)
	{
		forAll(groupSpecies[i],j)
		{
			if(speciesLambda.found(Y[groupSpecies[i][j]].name()))
			{
				lambdaDict = speciesLambda.subDict(Y[groupSpecies[i][j]].name());
				forAll(groupSpecies[i],j)
				{
					Lambda1[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
					Lambda2[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
					Lambda3[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
					Lambda4[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

					Lambda1[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda1")) * speciesLambdaUnit;
					Lambda2[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda2")) * speciesLambdaUnit;
					Lambda3[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda3")) * speciesLambdaUnit;
					Lambda4[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda4")) * speciesLambdaUnit;
				}
				break;
			}
			else 
			{
				if( j == (groupSpecies[i].size()-1) )
				{
					Info << "species conductivity log polynomial fitting coefficient for " << endl;
					Info << "(" << endl;
					forAll(groupSpecies[i],j)
					{
						Info << Y[groupSpecies[i][j]].name() << endl;
					}
					Info << ")" << endl;
					Info << " are not provided in constant/speciesLambda" <<	FatalErrorIn( "something wrong in your constant/speciesLambda") << abort(FatalError);
				}					
			}
		}
	}
}

Info << "species conductivity log polynomial fitting coefficient" << endl;
Info << "(" << endl;
forAll(Y,i)
{
		Info << Y[i].name() << tab << Lambda1[i].value() << tab << Lambda2[i].value() << tab << Lambda3[i].value() << tab << Lambda4[i].value() << endl;
}
Info << ")" << endl;
/*
forAll(Y,i)
{
  dictionary lambdaDict(speciesLambda.subDict(Y[i].name()));

  Lambda1[i].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
  Lambda2[i].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
  Lambda3[i].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
  Lambda4[i].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

  Lambda1[i] = readScalar(lambdaDict.lookup("Lambda1")) * specieslambdaUnit;
  Lambda2[i] = readScalar(lambdaDict.lookup("Lambda2")) * specieslambdaUnit;
  Lambda3[i] = readScalar(lambdaDict.lookup("Lambda3")) * specieslambdaUnit;
  Lambda4[i] = readScalar(lambdaDict.lookup("Lambda4")) * specieslambdaUnit;
}
*/