# Thioamide PolyProline PyRosetta
 PyRosetta Simulations of Thioamide FRET

## Algorithm Descriptions
### Thioamide_PolyProline.py
Simulates polyproline peptides containing fluorophore and backbone/sidechain thioamide pairs for comparisons to experimental FRET data.

### PolyProline_Data_Analysis.py
Performs analysis of outputs from Thioamide_PolyProline.py to compute fluorescence quenching via FRET, Dexter and distance-dependent quenching mechanisms.

### CaM_pOCNC_Thioamide.py
Simulates CaM_pOCNC complexes fluorophore and backbone thioamide pairs for comparisons to experimental FRET data.

### CaM_Data_Analysis.py
Performs analysis of outputs from CaM_pOCNC_Thioamide.py to compute fluorescence quenching via FRET, Dexter and distance-dependent quenching mechanisms.

### Compute_FRET_Interchanged.py
Performs all FRET calculations in the CaM_pOCNC_Thioamide.py script.

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


### Required Files
Since both backbone and sidechain thioamide simulations are run out of the same script all of the following files are required for runs:
```
TBL.params - Rosetta Residue Params file for parameterization of the backbone thioamide containing Leu amino acid
lys_thioacetyl.txt - Rosetta Residue Patch file for parameterizaiton of the sidechain, thioacetyl Lys amino acid
phe_cyanated.txt - Rosetta Residue Patch file for parameterization of the p-cyanophenylalanine fluorescent amino acid
thioamideN.txt - Rosetta Residue Patch file adjustment of the proceeding residue amide nitrogen charge to match that from Gaussian simulations of the thioamide
```
All files can be found in the params folder
## Thioamide PolyProline Simulations
### Running Thioamide_PolyProline
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
### Analyzing the PolyProline Output Data
Data analysis is performed with the PolyProline_Data_Analysis.py script as follows:
```
run PolyProline_Data_Analysis.py -F CNF -D Thioamide_Polyproline_BB -S RndChi_InterCIS_TermCIS_CHPiCIS
```
The flags are detailed below
```
-F --Donor_Fluor 3-Letter code for Donor fluorophore (CNF=cyanophenylalanine, TYR=tyrosine, TRP=tryptophan
-D --Directory_Name Name of the directory that houses outputs from all PolyProlines
-S --Sampling_Name  Name of the sampling scheme used. If no scheme was specified do not use this flag
```
Analysis script computes FRET using the PDA approximation, instantaneous kappa-squared and TrESP methods along with quenching due to distance-dependent quenching and Dexter along with all intermediate parameters.

## CaM/pOCNC Simulations
### Running CaM_pOCNC_Thioamide
The simulation can be run as follows, which represents the default simulation
```
run CaM_pOCNC_Thioamide.py -STRUCT_NUM 1 -RUB_NUM 10000 -OUT_NUM 100
```
Followed by the following flags detailed below
```
-STRUCT_NUM Identifies where the thioamide and fluorophore will be placed within the peptide/protein as detailed in the array on line 50.
-RUB_NUM Number of attempted Backrub/Sidechain Moves prior to structural output.
-OUT_NUM Number of output structures
-show Show simulation process in the command-line
-dump Dump the resultatnt structures to output PDBs
```
Note that the Compute_FRET_Interchanged.py script is required for running Cam/pOCNC simulations as this script specific performs the EFRET computes
### Analyzing the CaM/pOCNC Output Data
Data analysis is performed with the CaM_Data_Analysis.py script as follows:
```
run CaM_Data_Analysis.py -S Name_of_Cam_pOCNC_Thioamide_output_score_file.sc
```
Analysis script computes FRET using the PDA approximation, instantaneous kappa-squared and TrESP methods along with quenching due to distance-dependent quenching and Dexter along with all intermediate parameters.