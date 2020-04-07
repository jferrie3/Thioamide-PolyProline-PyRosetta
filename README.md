# Thioamide PolyProline PyRosetta
 PyRosetta Simulations of Thioamide FRET

## Algorithm Descriptions
### Thioamide_PolyProline.py


### Analyze_Thioamide_Polyproline.py


### Installation Guide
__Operating System:__ Linux (64-bit)

__Programming Langauge:__ Python
This code was specifically written and tested in Python3.6 (python3.6.8)
	
__Required Python Packages:__
- PyRosetta
	- This was specifically written and tested with PyRosetta 2019.17+release.2cb3f3a py36_0. Additional information can be found on the [PyRosetta](http://www.pyrosetta.org/) website. This specific PyRosetta build can be downloaded [here](http://www.pyrosetta.org/dow) after obtaining a [license](https://els.comotion.uw.edu/express_license_technologies/pyrosetta)
- biopython (1.73)
- numpy (1.14.5)
- scipy (1.1.0)

__Anaconda Environment:__
An anaconda environment containing all necessary packages can be found in the anaconda folder. Build time for this Anaconda environment takes on the order of mintues to hours depending on the number of processors used and is largely dependent on the PyRosetta build time. On a normal computer this is expected to take ~30 minutes. With this file you can generate a local version of this environment using the command:

```conda env create -f lion.yml```

## Simulation Times


## Required Files
Since both backbone and sidechain thioamide simulations are run out of the same script all of the following files are required for runs:
```
TBL.params - Rosetta Residue Params file for parameterization of the backbone thioamide containing Leu amino acid
lys_thioacetyl.txt - Rosetta Residue Patch file for parameterizaiton of the sidechain, thioacetyl Lys amino acid
phe_cyanated.txt - Rosetta Residue Patch file for parameterization of the p-cyanophenylalanine fluorescent amino acid
thioamideN.txt - Rosetta Residue Patch file adjustment of the proceeding residue amide nitrogen charge to match that from Gaussian simulations of the thioamide
```
All files can be found in the Required_Rosetta_Files folder
## Running Thioamide_PolyProline
The simulation can be run as follows, which represents the default simulation
```
run Thioamide_PolyProline.py -T BB -F CNF -PNUM 2 -Rnd_Chi False -Inter_CIS False -Term_CIS False -BB_NUM 100 -ROT_NUM 1000 -RUB_NUM 50 -show False -dump False
```
Followed by the following flags detailed below
```
-T --Simulation_Type Two letter code for switching between backbone (BB) and sidechain (SC) thioamide simulations
-F --Donor_Fluor 3-Letter code for Donor fluorophore (CNF=cyanophenylalanine, TYR=tyrosine, TRP=tryptophan
-PNUM --Number_of_Prolines Integer for number of interstitial prolines
-BB_NUM --Number_of_BB_Confs Number of tested backbone conformers
-ROT_NUM --Number_of_Rotamer_Confs Number of tested rotamers/side-conformers tested for each backbone geometry
-RUB_NUM --Number_of_Backrub_Moves Number of attempted Backrub/Sidechain Moves prior to structural output.
-show --Show_Results Number of tested rotamers/side-conformers tested for each backbone geometry
-dump --Dump_PDBs Dump resultant structures to output PDBs
-Term_CIS --Terminal_CIS False,help='Standard Internal 3% and Terminal 10% cis-Pro proabilities accounted for
-CHPI_CIS --CHPi_CIS Chromophore specific CH-Pi induced cisPro probabilities accounted for
-Inter_CIS --Internal_CIS Probability of Sampling Internal CIS residues at 2%
-Rnd_Chi --Random_Chi Effectively randomize the chi angles of terminal residues.
-CNF_TYR --CNF_as_TYR cis and trans probabilities from TYR rather than PHE
```
## Analyzing the Output Data
