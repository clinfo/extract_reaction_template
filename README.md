# extract_reaction_template
Extract reaction templates from chemical reaction data.

## Requirements
**Java: 1.8**  
**ChemAxon API: 18.18.0**  

## Installation
Refer to the link below.
* [**ChemAxon**](https://docs.chemaxon.com/display/docs/Installation+Guide)

## Usage
1. Donwload [USPTO reaction data](https://figshare.com/articles/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873) extracted by Lowe.
2. Get reaction smiles and remove duplicate reaction.
```bash
tail -n +2 1976_Sep2016_USPTOgrants_smiles.rsmi | awk -F '\t' '{print $1}' | sort | uniq > 1976_Sep2016_USPTOgrants_smiles_uniq.smi
```
3. Extract reaction templates.
```bash
javac ExtractReactionTemplate.java
java ExtractReactionTemplate -f smarts -i 1976_Sep2016_USPTOgrants_smiles_uniq.smi -o PATH_TO_OUTPUT
```
The output file contains; 
- product
- row reaction
- reaction template including the reaction center atoms
- reaction template including the reaction center atoms and 1-neighbor atoms

### Example
```bash
javac ExtractReactionTemplate.java
java ExtractReactionTemplate -f smarts -i sample.smi -o sample.txt
```


## Information
- [Reaction Mapping](https://docs.chemaxon.com/display/docs/Reaction+Mapping)
- [SMILES](https://docs.chemaxon.com/display/docs/SMILES)

## Reference
Shoichi Ishida , Kei Terayama,  Ryosuke Kojima, Kiyosei Takasu, Yasushi Okuno,
[Prediction and Interpretable Visualization of Synthetic Reactions Using Graph Convolutional Networks](http://dx.doi.org/10.26434/chemrxiv.8343995),
ChemRxiv. [DOI: 10.26434/chemrxiv.8343995]
