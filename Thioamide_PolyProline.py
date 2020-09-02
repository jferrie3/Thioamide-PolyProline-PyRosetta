# All probabilities used in the following simulation were selected from the following works
# According to to Dasgupta et.al. FEBS Letters 2007
## cisPro-cisPro-Phe -> 29% and Pro-cisPro-Phe -> 18%
## Phe chi1 -39.54 (gauche+) rotamer degrees for Pro(i)-Phe(i+2) interaction
## for cisPro-Phe, Phe chi1 g+ rotamer (N-CA-CB-CG -> 60 deg) (C-CA-CB-CG-120 -> -60)
# According to Ganguly et.al. JACS 2012
## cisPro-cisPro-Tyr (21.5%) and Pro-cisPro-Tyr (24.3%)
## prefered chi1 angle of -60
## aprox. 96(94) of cis-cis and 88(94) of trans-cis isomers showed CH-pi interaction
# According to Nardi et. al. J. Biomol. NMR 2000
## Terminal Tyr phi is -160/-80 in transPro and -155/-85 in cisPro
## psi is coil in both
## chi1 is g+(48%)/t(7%)/g-(45%) in transPro and g+(84%)/t(13%)/g-(3%) in cisPro
# According to Dubrack g+ -> 60 (-60), t -> 180 (60), g- -> -60 (180)
## By eye phi=-155, X1=60.0 X2=180.0 looks correct for C-H interaction
# Variables for quantum yield and J-overlaps are from: Cnf-> Goldberg et. al. JACS, Tyr-> Goldberg et. al. ChemComm, Trp-> Goldberg et. al. ChemComm
# Lifetimes from: Cnf-> Goldberg et. al., Tyr->Guzow et. al. J. Phys. Chem. B. 2004, Trp-> Ghisaidoobe, et. al. Int. J. Mol. Sci. 2014

#Setting Up Required Imports
from pyrosetta import *
init(extra_options =  '-relax::minimize_bond_angles True -packing:ex1:level 7 -packing:ex2:level 7 -packing:ex3:level 7 -packing:ex4:level 7 -extrachi_cutoff 0 -extra_res_fa ./TBL.params -extra_patch_fa ./lys_thioacetyl.txt -extra_patch_fa ./phe_cyanated.txt -extra_patch_fa ./thioamideN.txt')
from math import *
from random import random as rnd
from random import randint
from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue
from pyrosetta.rosetta.core.scoring import CA_rmsd
import numpy as np
import argparse

FOLDERNUMBER = 2
pmm = PyMOLMover()

# Argument Parsing
parser = argparse.ArgumentParser(description='Program')
parser.add_argument('-T', '--Simulation_Type', action='store', type=str, required=True, default='BB',
	help='Two letter code for switching between backbone (BB) and sidechain (SC) thioamide simulations')
parser.add_argument('-F', '--Donor_Fluor', action='store', type=str, required=True,
	help='3-Letter code for Donor fluorophore (CNF=cyanophenylalanine, TYR=tyrosine, TRP=tryptophan')
parser.add_argument('-PNUM', '--Number_of_Prolines', action='store', type=int, required=True,
	help='Integer for number of interstitial prolines')
parser.add_argument('-BB_NUM', '--Number_of_BB_Confs', action='store', type=int, required=False, default=100,
	help='Number of tested backbone conformers')
parser.add_argument('-ROT_NUM', '--Number_of_Rotamer_Confs', action='store', type=int, required=False, default = 1000,
	help='Number of tested rotamers/side-conformers tested for each backbone geometry')
parser.add_argument('-RUB_NUM', '--Number_of_Backrub_Moves', action='store', type=int, required=False, default = 50,
	help='Number of attempted Backrub/Sidechain Moves prior to structural output.')
parser.add_argument('-show', '--Show_Results', action='store', type=bool, required=False, default = False,
	help='Number of tested rotamers/side-conformers tested for each backbone geometry')
parser.add_argument('-dump', '--Dump_PDBs', action='store', type=bool, required=False, default = False,
	help='Dump resultant structures to output PDBs')
parser.add_argument('-Term_CIS', '--Terminal_CIS', action='store', type=bool, required=False, default = False,
	help='Standard Internal 3% and Terminal 10% cis-Pro proabilities accounted for')
parser.add_argument('-CHPI_CIS', '--CHPi_CIS', action='store', type=bool, required=False, default = False,
	help='Chromophore specific CH-Pi induced cisPro probabilities accounted for')
parser.add_argument('-Inter_CIS', '--Internal_CIS', action='store', type=bool, required=False, default = False,
	help='Probability of Sampling Internal CIS residues at 2%')
parser.add_argument('-Rnd_Chi', '--Random_Chi', action='store', type=bool, required=False, default = False,
	help='Effectively randomize the chi angles of terminal residues')
parser.add_argument('-CNF_TYR', '--CNF_as_TYR', action='store', type=bool, required=False, default = False,
	help='cis and trans probabilities from TYR rather than PHE')
args = parser.parse_args()

	
# Varables for Donor Fluorophore
donor_fluor_single_letter = ''
cis_cis_prob = 0.0; trans_cis_prob = 0.0; cis_trans_prob = 0.0; trans_trans_prob = 0.0
cc_chpi_prob = 0.0; tc_chpi_prob = 0.0
sc_cis_prob = 0.0; sc_trans_prob = 0.0

# Setting Specifics for each Donor
## cisPro probabilities from Ganguly et.al. JACS 2012
if args.Donor_Fluor == "CNF":
	donor_fluor_single_letter = 'F[PHE_p:cyanated]'
	cis_cis_prob = 0.186; trans_cis_prob = 0.171; cis_trans_prob = 0.095; trans_trans_prob = 0.548
	cc_chpi_prob = 0.94; tc_chpi_prob = 0.78
	sc_cis_prob = 0.23; sc_trans_prob = 0.77; sc_cisrot_prob = 0.62 # from Reimer et. al. JMB 1998 or 0.17 from Wu et. al. Biopolymers 1998
	if args.CNF_as_TYR == True:
		cis_cis_prob = 0.215; trans_cis_prob = 0.243; cis_trans_prob = 0.087; trans_trans_prob = 0.455
		cc_chpi_prob = 0.96; tc_chpi_prob = 0.88
		sc_cis_prob = 0.24; sc_trans_prob = 0.76; sc_cisrot_prob = 0.67 # or cis = 0.2
if args.Donor_Fluor == "TYR":
	donor_fluor_single_letter = 'Y'
	cis_cis_prob = 0.215; trans_cis_prob = 0.243; cis_trans_prob = 0.087; trans_trans_prob = 0.455
	cc_chpi_prob = 0.96; tc_chpi_prob = 0.88
	sc_cis_prob = 0.24; sc_trans_prob = 0.76; sc_cisrot_prob = 0.67 # or cis = 0.2
if args.Donor_Fluor == "TRP":
	donor_fluor_single_letter = 'W'
	cis_cis_prob = 0.285; trans_cis_prob = 0.166; cis_trans_prob = 0.065; trans_trans_prob = 0.484	
	cc_chpi_prob = 0.97; tc_chpi_prob = 0.81
	sc_cis_prob = 0.377; sc_trans_prob = 0.623; sc_cisrot_prob = 0.79 # or cis = 0.25

	
# Setting variable values for FRET calculations
QuantumYield = 0.0
JOverInt = 0.0
tau_d = 0.0 #seconds

if args.Donor_Fluor == "CNF":
	QuantumYield = 0.11
	JOverInt = 7*10**12
	tau_d = 7.0*10**(-9)
	if args.Simulation_Type == 'SC':
		QuantumYield = 0.11
		JOverInt = 2.37*10**12
		tau_d = 7.0*10**(-9)
if args.Donor_Fluor == "TYR":
	QuantumYield = 0.14
	JOverInt = 2.80*10**12
	tau_d = 3.4*10**(-9)
if args.Donor_Fluor == "TRP":
	QuantumYield = 0.13
	JOverInt = 1.63*10**9
	tau_d = 3.1*10**(-9)
	if args.Simulation_Type == 'SC':
		QuantumYield = 0.13
		JOverInt = 7.76*10**10
		tau_d = 3.1*10**(-9)

# Defining Custom Angle Adjustments for Experimental Matching
## Setting the structure to PPII
def reset_PPII(pose):
	for i in range (1, end+1):
		pose.set_phi(i, -78.0)
		pose.set_psi(i, 149)
		pose.set_omega(i, 180)	
	
## Probabilistic introduction of internal cisPro
def switch_internal(pose):
	cutoff_pro_num = end - 1
	start_pro_num = 1
	if args.Terminal_CIS == True:
		cutoff_pro_num = end - 2
		start_pro_num = 2
	if args.CHPi_CIS == True:
		cutoff_pro_num = end - 3
		if args.Simulation_Type == 'SC':
			cutoff_pro_num = end - 2
		start_pro_num = 2
	for i in range(start_pro_num, cutoff_pro_num):
		if 150 < pose.omega(i) < 210 or -210 < pose.omega(i) < -150:
			if rnd() < 0.02:
				pose.set_omega(i, 0)
				if rnd() <= 0.76:
					pose.set_phi(i, -76)
					pose.set_psi(i, 159)
				else:
					pose.set_phi(i, -86)
					pose.set_psi(i, 1)
		else:
			if rnd() < 0.98:
				pose.set_omega(i, 180)
				pose.set_phi(i, -78.0)
				pose.set_psi(i, 149)
			else:
				if rnd() <= 0.76:
					if pose.phi(i) < -81 or pose.phi(i) > -71:
						pose.set_phi(i, -76)
						pose.set_psi(i, 159)
				else:
					if pose.phi(i) < -91 or pose.phi(i) > -81:
						pose.set_phi(i, -86)
						pose.set_psi(i, 1)

##	Probabilistic introduction of N-terminal cisPro
def switch_cis_terminal(pose):
	if 150 < pose.omega(1) < 210 or -210 < pose.omega(1) < -150:
		cis_prob = 0.1 ## According to Reimer et. al. JMB 1998 is at max 0.12 for Leu-PP
		if rnd() < cis_prob:
			pose.set_omega(1, 0)
	else:
		trans_prob = 0.9
		if rnd() < trans_prob:
			pose.set_omega(1, 180)
	if 150 < pose.omega(end-2) < 210 or -210 < pose.omega(end-2) < -150:
		if rnd() < 0.1:
			pose.set_omega(end-2, 0)
	else:
		if rnd() < 0.9:
			pose.set_omega(end-2, 180)
			
## Probabilistic introduction of cisPro-Pro-Xxx, Pro-cisPro-Xxx and Pro-Pro-cisXxx
def switch_cis_CHPi(pose):
	if args.Simulation_Type == 'SC':
		## Yielding Aro-cisPro
		if 150 < pose.omega(1) < 210 or -210 < pose.omega(1) < -150: # trans
			if rnd() <= sc_cis_prob:
				pose.set_omega(1, 0)
				pose.set_phi(2, -76)
				pose.set_psi(2, 159)
		## Yielding Aro-transPro
		else: # cis
			if rnd() <= sc_trans_prob:
				pose.set_omega(1, 180)
				pose.set_phi(2, -78.0)
				pose.set_psi(2, 149)
	# Switching between transPro-Pro-Phe and cisPro-Pro-Phe
	else:
		## Yielding CisPro
		if 150 < pose.omega(end-3) < 210 or -210 < pose.omega(end-3) < -150: # trans
			if rnd() <= cis_cis_prob + cis_trans_prob:
				pose.set_omega(end-3, 0)
				if rnd() <= 0.74:
					pose.set_phi(end-3, -76)
					pose.set_psi(end-3, 159)
				else:
					pose.set_phi(end-3, -86)
					pose.set_psi(end-3, 1)
		## Yielding TransPro		
		else: #cis
			if rnd() <= trans_cis_prob + trans_trans_prob:
				pose.set_omega(end-3, 180)
				pose.set_phi(end-3, -78.0)
				pose.set_psi(end-3, 149)
			else:
				if rnd() <= 0.74:
					if pose.phi(end-3) < -81 or pose.phi(end-3) > -71:
						pose.set_phi(end-3, -76)
						pose.set_psi(end-3, 159)
				else:
					if pose.phi(end-3) < -91 or pose.phi(end-3) > -81:
						pose.set_phi(end-3, -86)
						pose.set_psi(end-3, 1)
		# Switching between Pro-transPro-Phe and Pro-cisPro-Phe
		if -30 < pose.omega(end-3) < 30 or 330 < pose.omega(end-3) < 390: #cis
			if 150 < pose.omega(end-2) < 210 or -210 < pose.omega(end-2) < -150: #trans
				if rnd() < cis_cis_prob/(cis_trans_prob + cis_cis_prob):
					pose.set_omega(end-2, 0)
					if rnd() <= 0.74:
						pose.set_phi(end-2, -76)
						pose.set_psi(end-2, 159)
					else:
						pose.set_phi(end-2, -86)
						pose.set_psi(end-2, 1)
			else: #cis 
				if rnd()< cis_trans_prob/(cis_trans_prob + cis_cis_prob):
					pose.set_omega(end-2, 180)
					pose.set_phi(end-2, -78.0)
					pose.set_psi(end-2, 149)
				else:
					if rnd() <= 0.74:
						if pose.phi(end-2) < -81 or pose.phi(end-2) > -71:
							pose.set_phi(end-2, -76)
							pose.set_psi(end-2, 159)
					else:
						if pose.phi(end-2) < -91 or pose.phi(end-2) > -81:
							pose.set_phi(end-2, -86)
							pose.set_psi(end-2, 1)
		else: #trans
			if 150 < pose.omega(end-2) < 210 or -210 < pose.omega(end-2) < -150: #trans
				if rnd() < trans_cis_prob/(trans_cis_prob + trans_trans_prob) :#0.18:
					pose.set_omega(end-2, 0)
					if rnd() <= 0.74:
						if pose.phi(end-2) < -81 or pose.phi(end-2) > -71:
							pose.set_phi(end-2, -76)
							pose.set_psi(end-2, 159)
					else:
						if pose.phi(end-2) < -91 or pose.phi(end-2) > -81:
							pose.set_phi(end-2, -86)
							pose.set_psi(end-2, 1)
			else: #cis
				if rnd() < trans_trans_prob/(trans_cis_prob + trans_trans_prob):#0.82:
					pose.set_omega(end-2, 180)
					pose.set_phi(end-2, -78.0)
					pose.set_psi(end-2, 149)
				else:	
					if rnd() <= 0.74:
						if pose.phi(end-2) < -81 or pose.phi(end-2) > -71:
							pose.set_phi(end-2, -76)
							pose.set_psi(end-2, 159)
					else:
						if pose.phi(end-2) < -91 or pose.phi(end-2) > -81:
							pose.set_phi(end-2, -86)
							pose.set_psi(end-2, 1)
						

## Probabilistic adjustment of chi angles for N-terminal Residue
def randomize_Lchi(pose):
	if args.Simulation_Type == 'SC':
		if -30 < pose.omega(end-2) < 30 or 330 < pose.omega(end-2) < 390:
			if rnd()<=0.5:
				pose.set_phi(int(end), -155)
				pose.set_psi(int(end), randint(-180, 180))
			else:
				pose.set_phi(int(end), -85)
				pose.set_psi(int(end), randint(-180, 180))
		else:
			if rnd() < 0.5:
				pose.set_phi(int(end), -160)
				pose.set_psi(int(end), randint(-180, 180))
			else:
				pose.set_phi(int(end), -80)
				pose.set_psi(int(end), randint(-180, 180))
		for i in range(1, pose.residue(int(end)).nchi()):
			pose.set_chi(i, int(end), randint(-180, 180))
		pose.set_chi(pose.residue(int(end)).nchi(),int(end), 180)	
	else:	
		if pose.omega(2) == 0:
			if rnd() < 0.5:
				pose.set_phi(1, -170)
			else:
				pose.set_phi(1, -70)
		else:
			if rnd() < 0.5:
				pose.set_phi(1, -165)
			else:
				pose.set_phi(1, -70)
		pose.set_psi(1, 159 + np.random.normal(0, 10))		
		for i in range(1, pose.residue(1).nchi()+1):
			pose.set_chi(i, 1, randint(-180, 180))

## Probabilistic adjustment of chi angles for C-terminal Residue
### transPro specific probabilistic adjustment
def transPro_rndchi(pose):
	rnd_4=rnd()
	if rnd_4 <= 0.2/(0.2 + 0.15 + 0.1):
		pose.set_chi(1, int(end), -60)
		pose.set_chi(2, int(end), randint(-180, 180))
	elif rnd_4 <= 0.35/(0.2 + 0.15 + 0.1):
		pose.set_chi(1, int(end), 60)
		pose.set_chi(2, int(end), randint(-180, 180))	
	else:
		pose.set_chi(1, int(end), 180)
		pose.set_chi(2, int(end), randint(-180, 180))
		
### cisPro specific probabilistic adjustment
def cisPro_rndchi(pose):
	rnd_4=rnd()
	if rnd_4 <= 0.1/(0.1 + 0.1 + 0.5):
		pose.set_chi(1, int(end), 60)
		if pose.omega(end-3) == 0:
			if rnd() <= cc_chpi_prob:
				if rnd() <0.5:
					pose.set_chi(2, int(end), 0)
				else:
					pose.set_chi(2, int(end), 180)
			else:
				pose.set_chi(2, int(end), randint(-180,180))
		else:
			if rnd() <= tc_chpi_prob:
				if rnd() <0.5:
					pose.set_chi(2, int(end), 0)
				else:
					pose.set_chi(2, int(end), 180)
			else:		
				pose.set_chi(2, int(end), randint(-180,180))
	elif rnd_4 <= 0.2/(0.1 + 0.1 + 0.5):
		pose.set_chi(1, int(end), 180)
		pose.set_chi(2, int(end), randint(-180, 180))	
	else:
		pose.set_chi(1, int(end), -60)
		if pose.omega(end-3) == 0:
			if rnd() <= cc_chpi_prob:
				if rnd() <0.5:
					pose.set_chi(2, int(end), -90)
				else:
					pose.set_chi(2, int(end), 90)
			else:
				pose.set_chi(2, int(end), randint(-180,180))
		else:
			if rnd() <= tc_chpi_prob:
				if rnd() <0.5:
					pose.set_chi(2, int(end), -90)
				else:
					pose.set_chi(2, int(end), 90)
			else:		
				pose.set_chi(2, int(end), randint(-180,180))

### cisPro specific probabilistic adjustment for N-terminal aromatic residue
def sc_cisPro_rndchi(pose):
	rnd_4=rnd()
	split_prob = (1-sc_cisrot_prob)/2
	if rnd_4 <= split_prob/(split_prob + sc_cisrot_prob + split_prob):
		pose.set_chi(1, 1, 60)		
		pose.set_chi(2, 1, randint(-180,180))
	elif rnd_4 <= sc_cisrot_prob/(split_prob + sc_cisrot_prob + split_prob):
		pose.set_chi(1, 1, 180)
		pose.set_chi(2, 1, 60)	
	else:
		pose.set_chi(1, 1, -60)
		pose.set_chi(2, 1, randint(-180,180))

### transPro specific probabilistic adjustment for N-terminal aromatic residue
def sc_transPro_rndchi(pose):
	rnd_4=rnd()
	if rnd_4 <= 0.2/(0.2 + 0.15 + 0.1):
		pose.set_chi(1, 1, -60)
		pose.set_chi(2, 1, randint(-180, 180))
	elif rnd_4 <= 0.35/(0.2 + 0.15 + 0.1):
		pose.set_chi(1, 1, 60)
		pose.set_chi(2, 1, randint(-180, 180))	
	else:
		pose.set_chi(1, 1, 180)
		pose.set_chi(2, 1, randint(-180, 180))
		
### Comibination of cis/transPro chi for adjustment
def randomize_Fchi(pose):
	if args.CHPi_CIS == True:
		if args.Simulation_Type == 'SC':
			if -30 < pose.omega(1) < 30 or 330 < pose.omega(1) < 390:
				sc_cisPro_rndchi(pose)
				if rnd() < 0.5:
					pose.set_phi(1, -170)
				else:
					pose.set_phi(1, -70)
			else:
				sc_transPro_rndchi(pose)
				if rnd() < 0.5:
					pose.set_phi(1, -165)
				else:
					pose.set_phi(1, -70)
			pose.set_psi(1, 159 + np.random.normal(0, 10))		
		else:
			if -30 < pose.omega(end-2) < 30 or 330 < pose.omega(end-2) < 390:
				cisPro_rndchi(pose)
				if rnd()<=0.5:
					pose.set_phi(int(end), -155)
					pose.set_psi(int(end), randint(-180, 180))
				else:
					pose.set_phi(int(end), -85)
					pose.set_psi(int(end), randint(-180, 180))
			else:
				transPro_rndchi(pose)
				if rnd() < 0.5:
					pose.set_phi(int(end), -160)
					pose.set_psi(int(end), randint(-180, 180))
				else:
					pose.set_phi(int(end), -80)
					pose.set_psi(int(end), randint(-180, 180))
	else:
		if args.Simulation_Type == 'SC':
			if -30 < pose.omega(1) < 30 or 330 < pose.omega(1) < 390:
				sc_cisPro_rndchi(pose)
				if rnd() < 0.5:
					pose.set_phi(1, -170)
				else:
					pose.set_phi(1, -70)
			else:
				sc_transPro_rndchi(pose)
				if rnd() < 0.5:
					pose.set_phi(1, -165)
				else:
					pose.set_phi(1, -70)
			pose.set_psi(1, 159 + np.random.normal(0, 10))	
			for i in range(1, 3):
				pose.set_chi(i, 1, randint(-180, 180))
		else:	
			if -30 < pose.omega(end-2) < 30 or 330 < pose.omega(end-2) < 390:
				if rnd()<=0.5:
					pose.set_phi(int(end), -155)
					pose.set_psi(int(end), randint(-180, 180))
				else:
					pose.set_phi(int(end), -85)
					pose.set_psi(int(end), randint(-180, 180))
			else:
				if rnd() < 0.5:
					pose.set_phi(1, -160)
					pose.set_psi(int(end), randint(-180, 180))
				else:
					pose.set_phi(1, -80)
					pose.set_psi(int(end), randint(-180, 180))
			for i in range(1, 3):
				pose.set_chi(i, int(end), randint(-180, 180))	
			
# SETTING THE PYROSETTA PARAMETERS
## STARTING STRUCTURE
p2_seq_str = "X[TBL]P[PRO_p:ThioamideN]"
if args.Simulation_Type == 'SC':
	p2_seq_str = donor_fluor_single_letter + "P"
for add_p in range(args.Number_of_Prolines-1):
	p2_seq_str += "P"
if args.Simulation_Type == 'SC':
	p2_seq_str += "K[LYS_p:thioacetylated]"
else:	
	p2_seq_str += donor_fluor_single_letter	
p2 = pose_from_sequence(p2_seq_str, "fa_standard")
end = p2.total_residue()
donor_residue = end
thioamide_residue = 1
if args.Simulation_Type == 'SC':
	donor_residue = 1
	thioamide_residue = end
add_variant_type_to_pose_residue(p2, pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT, end)
add_variant_type_to_pose_residue(p2, pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT, 1)
reset_PPII(p2)
startingp = Pose()
startingp.assign(p2)
p2_bb = Pose()

## SCORE FUNCTION
SF1 = create_score_function("ref2015_cart")
SF1.set_weight(pyrosetta.rosetta.core.scoring.dslf_fa13, 0)

## MONTE CARLO OBJECTS
mc1 = MonteCarlo(p2, SF1, 1.0)
mc2 = MonteCarlo(p2, SF1, 10.0) # originally 1.5
if args.Simulation_Type == 'SC':
	mc2 = MonteCarlo(p2, SF1, 10.0)
	
## MOVEMAP
movemap=MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)
movemap.set(pyrosetta.rosetta.core.id.DOF_Type.THETA, True)

## PACKER TASK
task = standard_packer_task(p2)
task.restrict_to_repacking()

## SIDECHAIN ROTAMER MOVER
sidechainMC = pyrosetta.rosetta.protocols.simple_moves.sidechain_moves.SidechainMover()
sidechainMC.set_task(task)
sidechainMC.set_change_chi_without_replacing_residue(True)
sidechainMC.set_prob_random_pert_current(0.0) # 0.8
sidechainMC.set_prob_uniform(0.55) # 0.2
sidechainMC.set_prob_withinrot(0.0) # 0.0

## MINIMZATION MOVER
minmover= pyrosetta.rosetta.protocols.minimization_packing.MinMover()
minmover.movemap(movemap)
minmover.score_function(SF1)
minmover.min_type('lbfgs_armijo_nonmonotone') #Call particular minimization type (dfpmin_armijo_nonmonotone) originally linmin
minmover.tolerance(0.0001)
minmover.cartesian(False)
minmover.omega(False)

## BACKRUB PROTOCOL
### PROTOCOL OPTIONS
prob_rotamer_only = 0.25
prob_rotamer_post_backrub = 0.25
prob_rotamer_pbr = 0.25
prob_backbone_swap = 0.35
backrub_mover = pyrosetta.rosetta.protocols.backrub.BackrubMover()

### BACKRUB PROTOCOL
def Backrub_Mover(p):
	mc2.reset(p)
	for i in range(args.Number_of_Backrub_Moves):
		if rnd() < prob_backbone_swap:
			randomize_Lchi(p)
			randomize_Fchi(p)
		if rnd() < prob_rotamer_only:
			sidechainMC.apply(p)
			p.set_phi(int(end), randint(-180, 180))
		else:
			mc1.reset(p)
			backrub_mover.apply(p)
			mc1.boltzmann(p)
			if rnd() < prob_rotamer_post_backrub:
				sidechainMC.apply(p)
				p.set_phi(int(end), randint(-180, 180))
				if rnd() < prob_rotamer_pbr:
					sidechainMC.apply(p)
					p.set_phi(int(end), randint(-180, 180))
		mc2.boltzmann(p)			
	return p

# Definitions for Vector Processing
## Definition for Vector Averaging
def vectoravg(vec1,vec2):
	avgvec = pyrosetta.rosetta.numeric.xyzVector_double_t(((vec1.x+vec2.x)/2.0),((vec1.y+vec2.y)/2.0),((vec1.z+vec2.z)/2.0))
	return avgvec
## Definition for Weighted Vector Averaging
def weightedvectoravg(weight1,vec1,weight2,vec2):
	weightavgvec = pyrosetta.rosetta.numeric.xyzVector_double_t((((weight1*vec1.x)+(weight2*vec2.x))/(weight1+weight2)),(((weight1*vec1.y)+(weight2*vec2.y))/(weight1+weight2)),(((weight1*vec1.z)+(weight2*vec2.z))/(weight1+weight2)))
	return weightavgvec

# Definitions for Counters
## Variables for Counters
cis_trans_counter = np.array([0, 0, 0, 0]) #[cc, tc, ct, tt]
running_TRESPEFRET_SUM = 0
running_TRESPEFRET_SUM_1 = 0
running_TRESPEFRET_SUM_2 = 0
running_EFRETKsq_SUM = 0
running_EFRETKsq_SUM_1 = 0
running_EFRETKsq_SUM_2 = 0
running_EFRETKsq_SUM_W = 0
running_EFRET_SUM = 0
running_EFRET_SUM_1 = 0
running_EFRET_SUM_2 = 0
running_EFRET_SUM_W = 0
running_TRESP_count = 0

## Cis/Trans-Pro_Cis/Trans-Pro_Xxx Counting Function
def cis_trans_counting(pose, cis_trans_counter):
	if -30 < pose.omega(end-3) < 30 or 330 < pose.omega(end-3) < 390:
		if -30 < pose.omega(end-2) < 30 or 330 < pose.omega(end-2) < 390:
			cis_trans_counter[0] += 1
		else:
			cis_trans_counter[2] += 1
	else:
		if -30 < pose.omega(end-2) < 30 or 330 < pose.omega(end-2) < 390:
			cis_trans_counter[1] += 1
		else:
			cis_trans_counter[3] += 1		
	return cis_trans_counter

# Definitions for Calculating EFRET
## Calculation of Average EFRET (Distance exclusive variable)
def compute_EFRET(pose_in):
	## Compute the R0
	R_0 = 0.211*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
	## Thioamide	
	if args.Simulation_Type == 'SC':
		ST = pose_in.residue(int(end)).xyz("ST")
		C = pose_in.residue(int(end)).xyz("CH")
		N = pose_in.residue(int(end)).xyz("NZ")
	else:
		ST = pose_in.residue(1).xyz("ST")
		C = pose_in.residue(1).xyz("C")
		N = pose_in.residue(2).xyz("N")
	CNmp = vectoravg(C,N)
	## Donor Fluor Variables
	if args.Donor_Fluor	== 'CNF': 
		## Cyanophenylalanine (from Goldberg et. al)
		if args.Simulation_Type == 'SC':
			CE1 = pose_in.residue(1).xyz("CE1")
			CE2 = pose_in.residue(1).xyz("CE2")
		else:
			CE1 = pose_in.residue(end).xyz("CE1")
			CE2 = pose_in.residue(end).xyz("CE2")
		CEmp = vectoravg(CE1,CE2)
		R = (CEmp - CNmp).norm()
		EFRET = 1/(1+(R/R_0)**6)
		return EFRET
	if args.Donor_Fluor	== 'TYR': 
		## Tyrosine (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		if args.Simulation_Type == 'SC':
			CE1 = pose_in.residue(1).xyz("CD1")
			CE2 = pose_in.residue(1).xyz("CD2")
		else:
			CE1 = pose_in.residue(end).xyz("CD1")
			CE2 = pose_in.residue(end).xyz("CD2")
		CEmp = vectoravg(CE1,CE2)
		R = (CEmp - CNmp).norm()
		EFRET = 1/(1+(R/R_0)**6)
		return EFRET
	if args.Donor_Fluor	== 'TRP': 
		## Tryptophan (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		### First Transition 4.8 M-1 x cm-1 x 10^(-3)
		if args.Simulation_Type == 'SC':
			NE = pose_in.residue(1).xyz("NE1")
			CE3 = pose_in.residue(1).xyz("CE3")
		else:
			NE = pose_in.residue(end).xyz("NE1")
			CE3 = pose_in.residue(end).xyz("CE3")
		CEmp_1 = vectoravg(NE,CE3)
		### Second Transition 0.8 M-1 x cm-1 x 10^(-3)
		if args.Simulation_Type == 'SC':
			CG = pose_in.residue(1).xyz("CG")
			CZ2 = pose_in.residue(1).xyz("CZ2")
		else:
			CG = pose_in.residue(end).xyz("CG")
			CZ2 = pose_in.residue(end).xyz("CZ2")
		CEmp_2 = vectoravg(CG,CZ2)
		### Combined Transition
		CEmp_W = weightedvectoravg(4.8, CEmp_1, 0.8, CEmp_2)
		R_1 = (CEmp_1 - CNmp).norm()
		R_2 = (CEmp_2 - CNmp).norm()
		R_W = (CEmp_W - CNmp).norm()
		EFRET_1 = 1/(1+(R_1/R_0)**6)
		EFRET_2 = 1/(1+(R_2/R_0)**6)
		EFRET_W = 1/(1+(R_W/R_0)**6)
		return EFRET_1, EFRET_2, EFRET_W

## Calculation of Dipole Corrected EFRET (Distance and k2 variable)
kappasquared = 0
### Compute KappaSquared from input donor and acceptor vectors
def compute_K2(pose_ST,pose_CNmp,pose_P1,pose_P2):
	### Computing Theta D
	V_pose_P1_pose_CNmp = pose_CNmp - pose_P1
	V_pose_P1_pose_p2 = pose_P2 - pose_P1
	Dot_pose_CNmp_pose_P1_pose_p2 = V_pose_P1_pose_CNmp.dot(V_pose_P1_pose_p2)
	l_pose_P1_pose_CNmp = V_pose_P1_pose_CNmp.norm()
	l_pose_P1_pose_p2 = V_pose_P1_pose_p2.norm()
	ThetaD = acos(Dot_pose_CNmp_pose_P1_pose_p2 / (l_pose_P1_pose_CNmp * l_pose_P1_pose_p2))
	### Computing Theta A
	V_pose_CNmp_ST = pose_ST - pose_CNmp
	V_pose_CNmp_pose_P1 = pose_P1 - pose_CNmp
	Dot_ST_pose_CNmp_pose_P1 = V_pose_CNmp_ST.dot(V_pose_CNmp_pose_P1)
	l_pose_CNmp_ST = V_pose_CNmp_ST.norm()
	l_pose_CNmp_pose_P1 = V_pose_CNmp_pose_P1.norm()
	ThetaA = acos(Dot_ST_pose_CNmp_pose_P1 / (l_pose_CNmp_ST * l_pose_CNmp_pose_P1))
	### Computing Theta DA	
	V_ST_pose_CNmp = pose_CNmp - pose_ST
	V_pose_CNmp_pose_P1 = pose_P1 - pose_CNmp
	V_pose_p2_pose_P1 = pose_P2 - pose_P1
	N1 = V_pose_CNmp_pose_P1.cross(V_ST_pose_CNmp) #cross vector of the two vectors forming a plane = normal vector of the plane
	N2 = V_pose_p2_pose_P1.cross(V_pose_CNmp_pose_P1)
	Dot_N1_N2 = N1.dot(N2)
	l_N1 = N1.norm()
	l_N2 = N2.norm()
	ThetaDA = acos(Dot_N1_N2 / (l_N1 * l_N2))
	K = (sin(ThetaD) * sin(ThetaA) * cos(ThetaDA)) - (2 * cos(ThetaD) * cos(ThetaA))
	K_sq = K ** 2
	return K_sq, ThetaD, ThetaA, ThetaDA
	
### Compute EFRET using KappaSquared from Dipoles
def compute_K2_EFRET(pose_in):
	## Compute the R0
	R_0 = 0.211*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
	## Thioamide	
	if args.Simulation_Type == 'SC':
		ST = pose_in.residue(int(end)).xyz("ST")
		C = pose_in.residue(int(end)).xyz("CH")
		N = pose_in.residue(int(end)).xyz("NZ")
	else:
		ST = pose_in.residue(1).xyz("ST")
		C = pose_in.residue(1).xyz("C")
		N = pose_in.residue(2).xyz("N")
	CNmp = vectoravg(C,N)
	## Donor Fluor Variables
	if args.Donor_Fluor	== 'CNF': 
		## Cyanophenylalanine (from Goldberg et. al)
		if args.Simulation_Type == 'SC':
			CE1 = pose_in.residue(1).xyz("CE1")
			CE2 = pose_in.residue(1).xyz("CE2")
		else:
			CE1 = pose_in.residue(end).xyz("CE1")
			CE2 = pose_in.residue(end).xyz("CE2")
		CEmp = vectoravg(CE1,CE2)
		R = (CEmp - CNmp).norm()
		K_sq, ThetaD, ThetaA, ThetaDA = compute_K2(ST,CNmp,CEmp,CE1)
		EFRETKsq = 1/(1+(R/(R_0*((K_sq/(2/3))**(1/6))))**6)
		return EFRETKsq, K_sq, ThetaD, ThetaA, ThetaDA, R
	if args.Donor_Fluor	== 'TYR': 
		## Tyrosine (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		if args.Simulation_Type == 'SC':
			CD1 = pose_in.residue(1).xyz("CD1")
			CD2 = pose_in.residue(1).xyz("CD2")
		else:
			CD1 = pose_in.residue(end).xyz("CD1")
			CD2 = pose_in.residue(end).xyz("CD2")
		CDmp = vectoravg(CD1,CD2)
		R = (CDmp - CNmp).norm()
		K_sq, ThetaD, ThetaA, ThetaDA = compute_K2(ST,CNmp,CD1,CD2)
		EFRETKsq = 1/(1+(R/(R_0*((K_sq/(2/3))**(1/6))))**6)
		return EFRETKsq, K_sq, ThetaD, ThetaA, ThetaDA, R
	if args.Donor_Fluor	== 'TRP': 
		## Tryptophan (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		### First Transition 4.8 M-1 x cm-1 x 10^(-3)
		if args.Simulation_Type == 'SC':
			NE = pose_in.residue(1).xyz("NE1")
			CE3 = pose_in.residue(1).xyz("CE3")
		else:
			NE = pose_in.residue(end).xyz("NE1")
			CE3 = pose_in.residue(end).xyz("CE3")
		CEmp_1 = vectoravg(NE,CE3)
		### Second Transition 0.8 M-1 x cm-1 x 10^(-3)
		if args.Simulation_Type == 'SC':
			CG = pose_in.residue(1).xyz("CG")
			CZ2 = pose_in.residue(1).xyz("CZ2")
		else:
			CG = pose_in.residue(end).xyz("CG")
			CZ2 = pose_in.residue(end).xyz("CZ2")
		CEmp_2 = vectoravg(CG,CZ2)
		### Combined Transition
		CEmp_W = weightedvectoravg(4.8,CEmp_1,0.8,CEmp_2)
		R_1 = (CEmp_1 - CNmp).norm()
		R_2 = (CEmp_2 - CNmp).norm()
		R_W = (CEmp_W - CNmp).norm()
		K_sq_1, ThetaD_1, ThetaA_1, ThetaDA_1 = compute_K2(ST,CNmp,NE,CE3)
		K_sq_2, ThetaD_2, ThetaA_2, ThetaDA_2 = compute_K2(ST,CNmp,CG,CZ2)
		K_sq_W = 4.8/(4.8+0.8)*K_sq_1 + 0.8/(4.8+0.8)*K_sq_2
		EFRETKsq_1 = 1/(1+(R_1/(R_0*((K_sq_1/(2/3))**(1/6))))**6)
		EFRETKsq_2 = 1/(1+(R_2/(R_0*((K_sq_2/(2/3))**(1/6))))**6)
		EFRETKsq_W = 1/(1+(R_W/(R_0*((K_sq_W/(2/3))**(1/6))))**6)
		return EFRETKsq_1, EFRETKsq_2, EFRETKsq_W, K_sq_1, ThetaD_1, ThetaA_1, ThetaDA_1, K_sq_2, ThetaD_2, ThetaA_2, ThetaDA_2, K_sq_W, R_1, R_2, R_W

## Calculation of TRESP EFRET (Excitonic Coupling)
## Constants
e_constant = 1.602*10**(-19) # Coulombs
D_Cm_Conv = 3.3356*10**(-30) # Debye to Coulomb meter conversion
def compute_TRESP_EFRET(pose_in):
	## Compute the R0
	R_0 = 0.211*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
	## All listed charges are in Atomic Units from Multiwfn fits from Gaussian CIS-aug-cc-PVTZ simulations
	## Thioamide Transition Charges
	if args.Simulation_Type == 'SC':
		ST_xyz = np.array(pose_in.residue(int(end)).xyz("ST"))
		C_xyz = np.array(pose_in.residue(int(end)).xyz("CH"))
		N_xyz = np.array(pose_in.residue(int(end)).xyz("NZ"))
	else:
		ST_xyz = np.array(pose_in.residue(1).xyz("ST"))
		C_xyz = np.array(pose_in.residue(1).xyz("C"))
		N_xyz = np.array(pose_in.residue(2).xyz("N"))
	#CA1_xyz = np.array(pose_in.residue(1).xyz("CA"))
	#CA2_xyz = np.array(pose_in.residue(2).xyz("CA"))
	atom_i_list_xyz = [ST_xyz, C_xyz, N_xyz]#, CA1_xyz, CA2_xyz]	
	ST_chg = 0.001093/sqrt(2)#0.012159/sqrt(2)
	C_chg = -0.001682/sqrt(2)#-0.006775/sqrt(2)
	N_chg = 0.000589/sqrt(2)#-0.010651/sqrt(2)
	#CA1_chg = -0.007636/sqrt(2)
	#CA2_chg = 0.002535/sqrt(2)
	atom_i_list_chg = [ST_chg, C_chg, N_chg]#, CA1_chg, CA2_chg]
	if args.Donor_Fluor	== 'CNF':
		## CNF Transition Charges from S0S1 of Cyanobenzene
		if args.Simulation_Type == 'SC':
			CE1_xyz = np.array(pose_in.residue(1).xyz("CE1"))
			CE2_xyz = np.array(pose_in.residue(1).xyz("CE2"))
			CD1_xyz = np.array(pose_in.residue(1).xyz("CD1"))
			CD2_xyz = np.array(pose_in.residue(1).xyz("CD2"))
			CG_xyz = np.array(pose_in.residue(1).xyz("CG"))
			CZ_xyz = np.array(pose_in.residue(1).xyz("CZ"))
			CT_xyz = np.array(pose_in.residue(1).xyz("CT"))
			NI_xyz = np.array(pose_in.residue(1).xyz("NI"))
		else:
			CE1_xyz = np.array(pose_in.residue(end).xyz("CE1"))
			CE2_xyz = np.array(pose_in.residue(end).xyz("CE2"))
			CD1_xyz = np.array(pose_in.residue(end).xyz("CD1"))
			CD2_xyz = np.array(pose_in.residue(end).xyz("CD2"))
			CG_xyz = np.array(pose_in.residue(end).xyz("CG"))
			CZ_xyz = np.array(pose_in.residue(end).xyz("CZ"))
			CT_xyz = np.array(pose_in.residue(end).xyz("CT"))
			NI_xyz = np.array(pose_in.residue(end).xyz("NI"))
		atom_j_list_xyz = [CE1_xyz, CE2_xyz, CD1_xyz, CD2_xyz, CG_xyz, CZ_xyz, CT_xyz, NI_xyz]
		CE1_chg = -0.119919/sqrt(2)
		CE2_chg = 0.119735/sqrt(2)
		CD1_chg = 0.019537/sqrt(2)
		CD2_chg = -0.018644/sqrt(2)
		CG_chg = -0.000169/sqrt(2)
		CZ_chg = -0.001698/sqrt(2)
		CT_chg = 0.001708/sqrt(2)
		NI_chg = -0.000549/sqrt(2)
		atom_j_list_chg = [CE1_chg, CE2_chg, CD1_chg, CD2_chg, CG_chg, CZ_chg, CT_chg, NI_chg]
		#### Compute J_TRESP
		J_TRESP = 0
		for atom_i in range(len(atom_i_list_xyz)):
			for atom_j in range(len(atom_j_list_xyz)):
				atom_ij_diff_xyz = np.linalg.norm(atom_j_list_xyz[atom_j]-atom_i_list_xyz[atom_i])*10**(-10)
				J_TRESP += (atom_i_list_chg[atom_i]*atom_j_list_chg[atom_j])/atom_ij_diff_xyz
		J_TRESP = (1)*J_TRESP*e_constant**2
		#### Compute Dipole Moments
		mu_A = 0
		for atom_i in range(len(atom_i_list_xyz)):
			mu_A += atom_i_list_chg[atom_i]*atom_i_list_xyz[atom_i]	
		mu_A = e_constant*np.linalg.norm(mu_A)*10**(-10) # in Coulomb meters
		mu_D = 0
		for atom_j in range(len(atom_j_list_xyz)):
			mu_D += atom_j_list_chg[atom_j]*atom_j_list_xyz[atom_j]
		mu_D = e_constant*np.linalg.norm(mu_D)*10**(-10) # in Coulomb meters
		TRESP_C = (3/2)*(((R_0*10**(-10))**6)/tau_d)*(1.33**(4))/((mu_A**2)*(mu_D**2))
		TRESP_k = J_TRESP**2*TRESP_C
		TRESPEFRET = 1/(1+(1/(tau_d*TRESP_k)))
		return TRESPEFRET
	if args.Donor_Fluor	== 'TYR':
		## TYR Transition Charges from S0S1 of Methylphenol	
		if args.Simulation_Type == 'SC':
			CE1_xyz = np.array(pose_in.residue(1).xyz("CE1"))
			CE2_xyz = np.array(pose_in.residue(1).xyz("CE2"))
			CD1_xyz = np.array(pose_in.residue(1).xyz("CD1"))
			CD2_xyz = np.array(pose_in.residue(1).xyz("CD2"))
			CG_xyz = np.array(pose_in.residue(1).xyz("CG"))
			CZ_xyz = np.array(pose_in.residue(1).xyz("CZ"))
			OH_xyz = np.array(pose_in.residue(1).xyz("OH"))
		else:
			CE1_xyz = np.array(pose_in.residue(end).xyz("CE1"))
			CE2_xyz = np.array(pose_in.residue(end).xyz("CE2"))
			CD1_xyz = np.array(pose_in.residue(end).xyz("CD1"))
			CD2_xyz = np.array(pose_in.residue(end).xyz("CD2"))
			CG_xyz = np.array(pose_in.residue(end).xyz("CG"))
			CZ_xyz = np.array(pose_in.residue(end).xyz("CZ"))
			OH_xyz = np.array(pose_in.residue(end).xyz("OH"))
		atom_j_list_xyz = [CE1_xyz, CE2_xyz, CD1_xyz, CD2_xyz, CG_xyz, CZ_xyz, OH_xyz]
		CE1_chg = 0.093182/sqrt(2)
		CE2_chg = -0.094158/sqrt(2)
		CD1_chg = 0.069516/sqrt(2)
		CD2_chg = -0.093729/sqrt(2)
		CG_chg = 0.007962/sqrt(2)
		CZ_chg = 0.011410/sqrt(2)
		OH_chg = 0.005818/sqrt(2)
		atom_j_list_chg = [CE1_chg, CE2_chg, CD1_chg, CD2_chg, CG_chg, CZ_chg, OH_chg]
		#### Compute J_TRESP
		J_TRESP = 0
		for atom_i in range(len(atom_i_list_xyz)):
			for atom_j in range(len(atom_j_list_xyz)):
				atom_ij_diff_xyz = np.linalg.norm(atom_j_list_xyz[atom_j]-atom_i_list_xyz[atom_i])*10**(-10)
				J_TRESP += (atom_i_list_chg[atom_i]*atom_j_list_chg[atom_j])/atom_ij_diff_xyz
		J_TRESP = (1)*J_TRESP*e_constant**2
		#### Compute Dipole Moments
		mu_A = 0
		for atom_i in range(len(atom_i_list_xyz)):
			mu_A += atom_i_list_chg[atom_i]*atom_i_list_xyz[atom_i]	
		mu_A = e_constant*np.linalg.norm(mu_A)*10**(-10) # in Coulomb meters
		mu_D = 0
		for atom_j in range(len(atom_j_list_xyz)):
			mu_D += atom_j_list_chg[atom_j]*atom_j_list_xyz[atom_j]
		mu_D = e_constant*np.linalg.norm(mu_D)*10**(-10) # in Coulomb meters
		TRESP_C = (3/2)*(((R_0*10**(-10))**6)/tau_d)*(1.33**(4))/((mu_A**2)*(mu_D**2))
		TRESP_k = J_TRESP**2*TRESP_C
		TRESPEFRET = 1/(1+(1/(tau_d*TRESP_k)))
		return TRESPEFRET
	if args.Donor_Fluor	== 'TRP':
		## TRP Transition Charges from S0S1 of Methyltryptophan
		if args.Simulation_Type == 'SC':
			CE2_xyz = np.array(pose_in.residue(1).xyz("CE2"))
			CE3_xyz = np.array(pose_in.residue(1).xyz("CE3"))
			CD1_xyz = np.array(pose_in.residue(1).xyz("CD1"))
			CD2_xyz = np.array(pose_in.residue(1).xyz("CD2"))
			CG_xyz = np.array(pose_in.residue(1).xyz("CG"))
			CZ2_xyz = np.array(pose_in.residue(1).xyz("CZ2"))
			CZ3_xyz = np.array(pose_in.residue(1).xyz("CZ3"))
			NE1_xyz = np.array(pose_in.residue(1).xyz("NE1"))
			CH2_xyz = np.array(pose_in.residue(1).xyz("CH2"))
		else:
			CE2_xyz = np.array(pose_in.residue(end).xyz("CE2"))
			CE3_xyz = np.array(pose_in.residue(end).xyz("CE3"))
			CD1_xyz = np.array(pose_in.residue(end).xyz("CD1"))
			CD2_xyz = np.array(pose_in.residue(end).xyz("CD2"))
			CG_xyz = np.array(pose_in.residue(end).xyz("CG"))
			CZ2_xyz = np.array(pose_in.residue(end).xyz("CZ2"))
			CZ3_xyz = np.array(pose_in.residue(end).xyz("CZ3"))
			NE1_xyz = np.array(pose_in.residue(end).xyz("NE1"))
			CH2_xyz = np.array(pose_in.residue(end).xyz("CH2"))
		atom_j_list_xyz = [CE2_xyz, CE3_xyz, CD1_xyz, CD2_xyz, CG_xyz, CZ2_xyz, CZ3_xyz, NE1_xyz, CH2_xyz]
		CE2_chg_1 = 0.001951/sqrt(2)
		CE3_chg_1 = -0.001083/sqrt(2) 
		CD1_chg_1 = 0.001676/sqrt(2)
		CD2_chg_1 = -0.000204/sqrt(2)
		CG_chg_1 = -0.000053/sqrt(2) 
		CZ2_chg_1 = 0.000242/sqrt(2) 
		CZ3_chg_1 = 0.000187/sqrt(2) 
		NE1_chg_1 = -0.002643/sqrt(2) 
		CH2_chg_1 = -0.000075/sqrt(2) 
		atom_j_list_chg_1 = [CE2_chg_1, CE3_chg_1, CD1_chg_1, CD2_chg_1, CG_chg_1, CZ2_chg_1, CZ3_chg_1, NE1_chg_1, CH2_chg_1]
		## TRP Transition Charges from S0S3 of Methyltryptophan
		CE2_chg_2 = -0.332447/sqrt(2)
		CE3_chg_2 = -0.498706/sqrt(2)
		CD1_chg_2 = 0.432104/sqrt(2)
		CD2_chg_2 = 0.375963/sqrt(2)
		CG_chg_2 = -0.215411/sqrt(2)
		CZ2_chg_2 = 0.525096/sqrt(2)
		CZ3_chg_2 = 0.413056/sqrt(2)
		NE1_chg_2 = -0.139179/sqrt(2)
		CH2_chg_2 = -0.560476/sqrt(2)
		atom_j_list_chg_2 = [CE2_chg_2, CE3_chg_2, CD1_chg_2, CD2_chg_2, CG_chg_2, CZ2_chg_2, CZ3_chg_2, NE1_chg_2, CH2_chg_2]
		#### Compute J_TRESP
		J_TRESP_1 = 0
		for atom_i in range(len(atom_i_list_xyz)):
			for atom_j in range(len(atom_j_list_xyz)):
				atom_ij_diff_xyz = np.linalg.norm(atom_j_list_xyz[atom_j]-atom_i_list_xyz[atom_i])*10**(-10)
				J_TRESP_1 += (atom_i_list_chg[atom_i]*atom_j_list_chg_1[atom_j])/atom_ij_diff_xyz
		J_TRESP_1 = (1)*J_TRESP_1*e_constant**2
		#### Compute Dipole Moments
		mu_A = 0
		for atom_i in range(len(atom_i_list_xyz)):
			mu_A += atom_i_list_chg[atom_i]*atom_i_list_xyz[atom_i]	
		mu_A = e_constant*np.linalg.norm(mu_A)*10**(-10) # in Coulomb meters
		mu_D_1 = 0
		for atom_j in range(len(atom_j_list_xyz)):
			mu_D_1 += atom_j_list_chg_1[atom_j]*atom_j_list_xyz[atom_j]
		mu_D_1 = e_constant*np.linalg.norm(mu_D_1)*10**(-10) # in Coulomb meters
		TRESP_C_1 = (3/2)*(((R_0*10**(-10))**6)/tau_d)*(1.33**(4))/((mu_A**2)*(mu_D_1**2))
		TRESP_k_1 = J_TRESP_1**2*TRESP_C_1
		TRESPEFRET_1 = 1/(1+(1/(tau_d*TRESP_k_1)))
		#### Compute J_TRESP
		J_TRESP_2 = 0
		for atom_i in range(len(atom_i_list_xyz)):
			for atom_j in range(len(atom_j_list_xyz)):
				atom_ij_diff_xyz = np.linalg.norm(atom_j_list_xyz[atom_j]-atom_i_list_xyz[atom_i])*10**(-10)
				J_TRESP_2 += (atom_i_list_chg[atom_i]*atom_j_list_chg_2[atom_j])/atom_ij_diff_xyz
		J_TRESP_2 = (1)*J_TRESP_2*e_constant**2
		#### Compute Dipole Moments
		mu_D_2 = 0
		for atom_j in range(len(atom_j_list_xyz)):
			mu_D_2 += atom_j_list_chg_2[atom_j]*atom_j_list_xyz[atom_j]
		mu_D_2 = e_constant*np.linalg.norm(mu_D_2)*10**(-10) # in Coulomb meters
		TRESP_C_2 = (3/2)*(((R_0*10**(-10))**6)/tau_d)*(1.33**(4))/((mu_A**2)*(mu_D_2**2))
		TRESP_k_2 = J_TRESP_2**2*TRESP_C_2
		TRESPEFRET_2 = 1/(1+(1/(tau_d*TRESP_k_2)))
		return TRESPEFRET_1,TRESPEFRET_2

# The Run Protocol
for bb_i in range (1, args.Number_of_BB_Confs+1):
	## Set up the backbone geometry
	p2.assign(startingp)
	reset_PPII(p2)
	if args.Terminal_CIS == True:
		switch_cis_terminal(p2)
	if args.CHPi_CIS == True:	
		switch_cis_CHPi(p2)
	if args.Internal_CIS == True:
		switch_internal(p2)
	if args.Random_Chi == True:	
		randomize_Lchi(p2)
		randomize_Fchi(p2)
	minmover.apply(p2)
	mc1.reset(p2)
	mc2.reset(p2)
	p2_bb.assign(p2)
	for rot_k in range(1, args.Number_of_Rotamer_Confs+1):
		## Sample Rotamers
		p2.assign(p2_bb)
		p2 = Backrub_Mover(p2)
		## Perform EFRET Computes and Dumps	
		running_TRESP_count += 1
		### For CNF and Tyr
		if args.Donor_Fluor == 'CNF' or args.Donor_Fluor == 'TYR':
			EFRET = compute_EFRET(p2)
			running_EFRET_SUM += EFRET
			running_EFRET_AVG = running_EFRET_SUM/running_TRESP_count
			EFRETKsq, K_sq, ThetaD, ThetaA, ThetaDA, R = compute_K2_EFRET(p2)
			running_EFRETKsq_SUM += EFRETKsq
			running_EFRETKsq_AVG = running_EFRETKsq_SUM/running_TRESP_count
			TRESPEFRET = compute_TRESP_EFRET(p2)
			running_TRESPEFRET_SUM += TRESPEFRET
			running_TRESP_AVG = running_TRESPEFRET_SUM/running_TRESP_count
			outf_name = args.Donor_Fluor + "_P" + str(args.Number_of_Prolines)
			if args.Random_Chi == True:
				outf_name += "_RndChi"
			if args.Internal_CIS == True:
				outf_name += "_InterCIS"
			if args.Terminal_CIS == True:
				outf_name += "_TermCIS"
			if args.CHPi_CIS == True:
				outf_name += "_CHPiCIS"
			outf_chi_name = outf_name + ".chis"	
			outf_name += ".score"
			outf = open(outf_name, 'a')
			outf_chi = open(outf_chi_name, 'a')
			pdb_out = args.Donor_Fluor + "_P" + str(args.Number_of_Prolines) + "_" + str(bb_i) + '_' + str(rot_k) + ".pdb"
			outf.write("%s\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (pdb_out, SF1(p2), CA_rmsd(startingp, p2), R, K_sq, ThetaD, ThetaA, ThetaDA, EFRET, EFRETKsq, TRESPEFRET, running_EFRET_AVG, running_EFRETKsq_AVG, running_TRESP_AVG))
			outf_chi.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (pdb_out, p2.residue(1).chi(1), p2.residue(1).chi(2), p2.residue(end).chi(1), p2.residue(end).chi(2), p2.residue(end).chi(3), p2.residue(end).chi(4), p2.residue(end).chi(5)))
			if args.Dump_PDBs and args.Dump_PDBs == True:
				p2.dump_pdb(pdb_out)
			outf.close()
			outf_chi.close()
			if args.Show_Results and args.Show_Results == True:
				pmm.apply(p2)
				print('Instantaneous Energy:' + str(SF1(p2)))
				print('Average EFRET:' + str(running_EFRET_AVG))
				print('Average EFRETKsq:' + str(running_EFRETKsq_AVG))
				print('Average TRESP EFRET:' + str(running_TRESP_AVG))
				print('[cis-cis, trans-cis, cis-trans, trans-trans]' + str(cis_trans_counting(p2, cis_trans_counter)/running_TRESP_count))
		### For Trp		
		if args.Donor_Fluor == 'TRP':
			EFRET_1, EFRET_2, EFRET_W = compute_EFRET(p2)
			running_EFRET_SUM_1 += EFRET_1
			running_EFRET_AVG_1 = running_EFRET_SUM_1/running_TRESP_count
			running_EFRET_SUM_2 += EFRET_2
			running_EFRET_AVG_2 = running_EFRET_SUM_2/running_TRESP_count
			running_EFRET_SUM_W += EFRET_W
			running_EFRET_AVG_W = running_EFRET_SUM_W/running_TRESP_count
			EFRETKsq_1, EFRETKsq_2, EFRETKsq_W, K_sq_1, ThetaD_1, ThetaA_1, ThetaDA_1, K_sq_2, ThetaD_2, ThetaA_2, ThetaDA_2, K_sq_W, R_1, R_2, R_W = compute_K2_EFRET(p2)
			running_EFRETKsq_SUM_1 += EFRETKsq_1
			running_EFRETKsq_AVG_1 = running_EFRETKsq_SUM_1/running_TRESP_count
			running_EFRETKsq_SUM_2 += EFRETKsq_2
			running_EFRETKsq_AVG_2 = running_EFRETKsq_SUM_2/running_TRESP_count
			running_EFRETKsq_SUM_W += EFRETKsq_W
			running_EFRETKsq_AVG_W = running_EFRETKsq_SUM_W/running_TRESP_count
			TRESPEFRET_1,TRESPEFRET_2 = compute_TRESP_EFRET(p2)
			running_TRESPEFRET_SUM_1 += TRESPEFRET_1
			running_TRESPEFRET_SUM_2 += TRESPEFRET_2
			running_TRESP_AVG_1 = running_TRESPEFRET_SUM_1/running_TRESP_count
			running_TRESP_AVG_2 = running_TRESPEFRET_SUM_2/running_TRESP_count
			outf_name = args.Donor_Fluor + "_P" + str(args.Number_of_Prolines)
			if args.Random_Chi == True:
				outf_name += "_RndChi"
			if args.Internal_CIS == True:
				outf_name += "_InterCIS"
			if args.Terminal_CIS == True:
				outf_name += "_TermCIS"
			if args.CHPi_CIS == True:
				outf_name += "_CHPiCIS"
			outf_name += ".score"
			outf = open(outf_name, 'a')
			pdb_out = args.Donor_Fluor + "_P" + str(args.Number_of_Prolines) + "_" + str(bb_i) + '_' + str(rot_k) + ".pdb"
			outf.write("%s\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (pdb_out, SF1(p2), CA_rmsd(startingp, p2), R_1, R_2, R_W, K_sq_1, ThetaD_1, ThetaA_1, ThetaDA_1, K_sq_2, ThetaD_2, ThetaA_2, ThetaDA_2, K_sq_W, EFRET_1, EFRET_2, EFRET_W, EFRETKsq_1, EFRETKsq_2, EFRETKsq_W, TRESPEFRET_1, TRESPEFRET_2, running_EFRET_AVG_1, running_EFRET_AVG_2, running_EFRET_AVG_W, running_EFRETKsq_AVG_1, running_EFRETKsq_AVG_2, running_EFRETKsq_AVG_W, running_TRESP_AVG_1, running_TRESP_AVG_2))
			if args.Dump_PDBs and args.Dump_PDBs == True:
				p2.dump_pdb(pdb_out)
			outf.close()
			if args.Show_Results and args.Show_Results == True:
				pmm.apply(p2)
				print('Instantaneous Energy:' + str(SF1(p2)))
				print('Average EFRET:' + str(running_EFRET_AVG_1) + ' ' + str(running_EFRET_AVG_2) + ' ' + str(running_EFRET_AVG_W))
				print('Average EFRETKsq:' + str(running_EFRETKsq_AVG_1) + ' ' + str(running_EFRETKsq_AVG_2) + ' ' + str(running_EFRETKsq_AVG_W))
				print('Average TRESP EFRET:' + str(running_TRESP_AVG_1) + ' ' + str(running_TRESP_AVG_2))
				print('[cis-cis, trans-cis, cis-trans, trans-trans]' + str(cis_trans_counting(p2, cis_trans_counter)/running_TRESP_count))

		