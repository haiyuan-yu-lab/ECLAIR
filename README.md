## Input Format
The pipeline only requires a tab-delimited text file with the pairs of interacting proteins that will be used to predict interfaces.
```
Q6VAB6	Q02750
Q6VAB6	P28482
Q6VAB6	P15056
Q6VAB6	P54646
```

## Output Format
Once the pipeline is done running, a CSV file will be generated per interaction in the chosen output folder.

| Column | Contents |
|--------|---------|
| P1 |	Protein Partner 1 (Alphabetical) |
| P2 | Protein Partner 2 (Alphabetical) |
| Prot |	Indicated whether this row refers to the first (0) or the second (1) protein partner |
| Pos	| 1-indexed Amino Acid Position that this row refers to |
| Res	| Expected Amino Acid at this position |
| TopClf |	Which of the 8 ECLAIR classifiers was used for this prediction |
| Pred	| Probability of this being an interface residue |
| Interface Prediction	| Binned probability of this being an interface residue (Very Low / Low / Medium / High / Very High) |


```
P1,P2,Prot,Pos,Res,TopClf,Pred,Interface Prediction
P15056,Q6VAB6,0,1,M,0,0.0,Very Low
P15056,Q6VAB6,0,2,A,0,0.0,Very Low
P15056,Q6VAB6,0,3,A,0,0.0,Very Low
P15056,Q6VAB6,0,4,L,0,0.0,Very Low
P15056,Q6VAB6,0,5,S,2,0.00169406074475,Very Low
P15056,Q6VAB6,0,6,G,2,0.00286604835082,Very Low
P15056,Q6VAB6,0,7,G,2,0.00219393392881,Very Low
P15056,Q6VAB6,0,8,G,2,0.000384057971014,Very Low
P15056,Q6VAB6,0,9,G,2,0.0,Very Low
```

## Running the script
The script requires a specific anaconda environment to run: ecalir.txt

Note that running the script is extremely resource intensive. As a rule of thumb, allocate 10GB of RAM per core used.


Code
The code below outlines the main script that orchestrates the data parsing subroutines for ECLAIR. As you can see in the code, the process is broken into 8 steps, each of which will run many different helper scripts.

The individual helper scripts can be used independently to generate individual features.
ECLAIR Pipeline Steps
* STEP 0	Code Setup	Sets up logging behavior which can be accessed at /home/jc2375/outputs/get_features.log
* STEP 1	List new UniProt Ids	Downloads the latest information for the UniProt IDs in the file and appends them to a cached uniprot_info.txt file
* STEP 2	ExPasy and Pfam features	Extracts information about biophysical features, isoforms, and Pfam domains
* STEP 3	Conservation Features	Generates multiple sequence alignment files and calculates both Jensen-Shannon Conservation and Coevolution via DCA and SCA
* STEP 4	MODBASE Models	Downloads existing MODBASE models and calculates surface residues
* STEP 5	PDB Structures	Downloads existing PDB models and calculates surface residues
* STEP 6	DOCKING	Runs several trials of molecular docking to generate paired structural features when available
* STEP 7	Make Interaction Files	Compiles all information parsed above into ECLAIR-compatible pandas dataframes
* STEP 8	Predict	Runs the ECLAIR classifier and saves results as either pickled dataframes or csv files
