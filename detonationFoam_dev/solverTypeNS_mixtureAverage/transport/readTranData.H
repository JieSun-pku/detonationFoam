//***************** S. Zhong 20201008******************//
// Considering that plenty of species have the same transport properties, therefore the transdata from one species can be used for other similar species//
// it can also save some computational time.
/*
e.g. 
NC7H16              2   459.6   6.253   0.0 0.0 1.0 !   TCPC                  
C7H15-1             2   459.6   6.253   0.0 0.0 1.0 !   WJP                   
C7H15                2   459.6   6.253   0.0 0.0 1.0 !   WJP       
*/
//First read the tran.data
IFstream speciesInfo(mesh.time().constant()/"trandat");
if(!speciesInfo.good())
{
	Info << FatalErrorIn( "cannot open trandat in constant") << abort(FatalError);
}
Foam::string   headerStr;
Foam::word    thisSpecies;
boolList foundSpecies(Y.size(), false);
//Foam::string   thisSpecies;
List<scalarList> trandata(Y.size());
forAll(trandata,i)
{
	trandata[i].setSize(6, 0.0);
}

speciesInfo.Foam::ISstream::getLine(headerStr);
while( !speciesInfo.eof() )
{
	speciesInfo >> thisSpecies;
	if(!speciesInfo.bad() )
	{
		for(label ii=0; ii< Y.size() ; ii++)
		{
			if(thisSpecies == Y[ii].name())
			{ 
				foundSpecies[ii] = true;
				speciesInfo >> trandata[ii][0] >> trandata[ii][1]  >> trandata[ii][2] >> trandata[ii][3]  >> trandata[ii][4]  >> trandata[ii][5]  ;
				Info <<"read the trandata of "<< thisSpecies  << tab<<  trandata[ii][0] << tab << trandata[ii][1]<< tab << trandata[ii][2]<< tab << trandata[ii][3]<< tab<<trandata[ii][4]<<  tab<<trandata[ii][5]<< tab<< endl;
			}
		}
	}
	speciesInfo.Foam::ISstream::getLine(headerStr);	
}

forAll(Y,i)
{
	if(!foundSpecies[i])
	{
		Info << " species trandata for " << Y[i].name() << " is not found in constant/trandat "<<	FatalErrorIn( "something wrong in your constant/trandat") << abort(FatalError);	
	}		
}
// group the species that have the same trandata
List<scalarList> groupSpecies(Y.size());
forAll(groupSpecies,i)
{
	groupSpecies[i].setSize(1,i);
}

forAll(trandata,i)
{
	for(label j = 0; j < i;j++)
	{
		if(i != j)
		{
			if(trandata[i] == trandata[j])
			{
				label currentSize =  groupSpecies[j].size();
				groupSpecies[j].resize(currentSize+1);
				groupSpecies[j][currentSize] = i;
				groupSpecies[i][0] = -1;
				break;
			}	
		}
	}
}

forAll(groupSpecies,i)
{
	if(groupSpecies[i].size() > 1)
	{
		Info << "trandata group species:  " << endl;
		Info << "(" << endl;
		forAll(groupSpecies[i],j)
		{
				Info << Y[groupSpecies[i][j]].name() << endl;
		}
		Info << ")" << endl;
	}
}

Info << "trandata reduced species list:  " << endl;
Info << "(" << endl;
forAll(groupSpecies,i)
{
	if(groupSpecies[i][0] != -1)
	{
		Info << Y[groupSpecies[i][0]].name() << endl;		
	}
}
Info << ")" << endl;
