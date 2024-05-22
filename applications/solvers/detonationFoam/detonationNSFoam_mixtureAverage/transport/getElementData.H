/* label  nSpecies = composition.Y().size();
List<List<specieElement>> specieComp(Y.size());
speciesTable speciesTab = thermo.composition().species();
const HashTable<List<specieElement>>& specComp = dynamic_cast<const multiComponentMixture<sutherlandTransport<species::thermo<janafThermo<perfectGas<specie>>,sensibleEnthalpy>>>&>(thermo).specieComposition();
forAll(specieComp,i)
{
    specieComp[i] = specComp[Y[i].member()];
}
wordList elementNames(4);

{
    elementNames[0] = "H";
    elementNames[1] = "C";
    elementNames[2] = "O";
    elementNames[3] = "N";
}

List<scalarField> spInfo(nSpecies);
forAll(spInfo,i)  spInfo[i].setSize(8,0);

forAll(Y, i)
{
    Info << "specie = " << Y[i].name() << "\t";
    forAll(elementNames, en)
    {
            forAll(specieComp[i], kk)
            {
                if(specieComp[i][kk].name() == elementNames[en])
                {
                    Info << specieComp[i][kk].name() << ", nAtoms = " << specieComp[i][kk].nAtoms() << "\t";
                    spInfo[i][en] = specieComp[i][kk].nAtoms();
                }
            }

    }

    Info << composition.Wi(i) << endl;
    spInfo[i][4] = composition.Wi(i);
}
 */