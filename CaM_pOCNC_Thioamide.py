## EXTRA IMPORT FOR FRET CALCULATION
from Compute_FRET_Interchanged import *

## STARTING UP PYROSETTA ##
from pyrosetta import * 
init(extra_options =  '-extrachi_cutoff 0 -extra_res_fa ./TBL.params -extra_res_fa ./TBF.params -extra_res_fa ./CNF.params -extra_patch_fa ./thioamideN.txt')

## IMPORTING FUNCTIONS ##
import csv
import sys
import argparse
from math import exp, log, pi, sqrt
from random import random as rnd
from random import randint
import numpy as np
import sys
from Bio.PDB import *
import argparse

## EXTRA ROSETTA IMPORTS ##
from pyrosetta.rosetta.protocols import backrub
from pyrosetta.rosetta.protocols.simple_moves import sidechain_moves
from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue
from pyrosetta.rosetta.core.scoring import CA_rmsd

pmm = PyMOLMover()

# Argument Parsing
parser = argparse.ArgumentParser(description='Program')
parser.add_argument('-STRUCT_NUM', '--Number_of_Input_Structure', action='store', type=int, required=True,
	help='Conformer number for Input PDB')
parser.add_argument('-RUB_NUM', '--Number_of_Backrub_Moves', action='store', type=int, required=True,
	help='Number of attempted Backrub/Sidechain Moves prior to structural output. Default = 10000')
parser.add_argument('-OUT_NUM', '--Number_of_Output_Structures', action='store', type=int, required=True,
	help='Number of output structures. Default = 100')
parser.add_argument('-show', '--Show_Results', action='store', type=bool, required=False,
	help='Show simulation progress in the command-line')
parser.add_argument('-dump', '--Dump_PDBs', action='store', type=bool, required=False,
	help='Dump resultant structures to output PDBs')
args = parser.parse_args()

## NEED TO SIMULATE
### CNF CaM/Peptide (QuantumYield/EQ)
#F'1 = res 153
# F*13/F'1 (0.003/0.11), F*13/L'11 (0.003/0.37), F*93/F'1 (0.014/0.79), F*93/L'11 (0.014/0.43), F*100/F'1(0.110/0.25), F*100/L'11(0.110/0.17), F*100/F'16 172 (0.110/0.10)
### Tyr CaM/Peptide (QuantumYield/EQ)
# Y139/F'1 (0.061/0.50), W100/F'1 (0.15/0.27), Y100/F'1(0.105/0.22), W139/F'1(0.29/0.58)

simulation_list = [['CNF', 12, 'TBF', 153, 0.11, 0.11],
	['CNF', 12, 'TBL', 163, 0.11, 0.37],
	['CNF', 92, 'TBF', 153, 0.014, 0.79],
	['CNF', 92, 'TBL', 163, 0.014, 0.43],
	['CNF', 99, 'TBF', 153, 0.110, 0.25],
	['CNF', 99, 'TBL', 163, 0.110, 0.17],
	['CNF', 99, 'TBF', 168, 0.110, 0.10],
	['TYR', 138, 'TBF', 153, 0.061, 0.50],
	['TRP', 99, 'TBF', 153, 0.15, 0.27],
	['TYR', 99, 'TBF', 153, 0.105, 0.22],
['TRP', 138, 'TBF', 153, 0.29, 0.58]]

## POSE AND MUTATION ##
"""How residues are numbered in PyMOL is different from how it is numbered in Rosetta.
Use $print pose.residue(5)$ to make sure that the residue numbers are correct."""
starting_pose = pose_from_pdb('Relaxed_1SY9_' + str(args.Number_of_Input_Structure) + '.pdb')
p = Pose()
p.assign(starting_pose)
mutant_start_pose = Pose()

def Mutate(pose, Fluor, FluorRes, ThioAA, ThioRes):
	AllPhe0 = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(99, "PHE")
	AllPhe1 = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(138, "PHE")
	AllPhe2 = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(166, "PHE")
	AllPhe3 = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(168, "PHE")
	MutThio = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(ThioRes, ThioAA)
	MutFluor = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(FluorRes, Fluor)
	AllPhe0.apply(pose) ; AllPhe1.apply(pose) ; AllPhe2.apply(pose)
	AllPhe3.apply(pose) ; MutThio.apply(pose) ; MutFluor.apply(pose)
	add_variant_type_to_pose_residue(pose, pyrosetta.rosetta.core.chemical.VariantType.SIDECHAIN_CONJUGATION, ThioRes+1)
			
## SCORE FUNCTION
SF1 = create_score_function("ref2015")
SF1.set_weight(pyrosetta.rosetta.core.scoring.dslf_fa13, 0)

## MONTE CARLO OBJECTS
mc1 = MonteCarlo(p, SF1, 1.0)
mc2 = MonteCarlo(p, SF1, 10.0)

## SET UP Movemap ##
movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)
#movemap.set(pyrosetta.rosetta.core.id.DOF_Type.THETA, True)


## PACKER TASK
# Residues 149, 150, 151, 152 are Calcium
task = standard_packer_task(starting_pose)
task.restrict_to_repacking()
'''vector = rosetta.utility.vector1_bool()
for res_i in range(1, p.total_residue()+1):
	if res_i == 149 or res_i == 150 or res_i == 151 or res_i == 152:
		vector.append(False)
	else:
		vector.append(True)
task.restrict_to_residues(vector)'''

## SIDECHAIN ROTAMER MOVER
sidechainMC = pyrosetta.rosetta.protocols.simple_moves.sidechain_moves.SidechainMover()
sidechainMC.set_task(task)
sidechainMC.set_change_chi_without_replacing_residue(True)
#sidechainMC.set_prob_random_pert_current(0.0) # 0.8
#sidechainMC.set_prob_uniform(0.55) # 0.2
#sidechainMC.set_prob_withinrot(0.0) # 0.0
packed_residue_list = sidechainMC.packed_residues()

## PACKER FOR MUTATIONS
pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(SF1, task)

## MINIMZATION MOVER
minmover= pyrosetta.rosetta.protocols.minimization_packing.MinMover()
minmover.movemap(movemap)
minmover.score_function(SF1)
minmover.min_type('lbfgs_armijo_nonmonotone') #Call particular minimization type (dfpmin_armijo_nonmonotone) originally linmin
minmover.tolerance(0.0001)
#minmover.cartesian(True)
minmover.omega(False)


## BACKRUB PROTOCOL
### PROTOCOL OPTIONS
prob_rotamer_only = 0.25
prob_rotamer_post_backrub = 0.25
prob_rotamer_pbr = 0.25
backrub_mover = pyrosetta.rosetta.protocols.backrub.BackrubMover()
backrub_mover.apply(p) ### backrub mover appplied to get the segment IDs

### SEGMENT IDs
residue_1 = []
residue_2 = []
backrub_seg_ids = []

for res1 in range(1, p.total_residue()+1, 1):
	for res2 in range(1, p.total_residue()+1, 1):
		residue_1.append(res1)
		residue_2.append(res2)
		res1ca = AtomID(2,res1)
		res2ca = AtomID(2,res2)
		backrub_seg_ids.append(backrub_mover.segment_id(res1ca,res2ca))
active_backrub_seg_ids = [i for i, seg_id in enumerate(backrub_seg_ids) if seg_id != 0]

#### DEFINED MOVER COMPONENTS ####
def backrub_sidechain_mover(p):
	last_backrub_seg_id = backrub_seg_ids.index(backrub_mover.last_segment_id())
	side_move_start_res = residue_1[last_backrub_seg_id]
	side_move_end_res = residue_2[last_backrub_seg_id]
	sidechain_movemap = MoveMap()
	is_pack_res = False
	res_attempt = 0
	sc_res_num = 0
	while is_pack_res == False:
		if res_attempt < 10:
			res_attempt = res_attempt + 1
			#print(res_attempt)
			sc_res_num = randint(side_move_start_res,side_move_end_res)
			if sc_res_num in packed_residue_list:
				is_pack_res = True
		else:
			sc_res_num_rand = randint(1,len(packed_residue_list))
			sc_res_num = packed_residue_list[sc_res_num_rand]	
			is_pack_res = True
	sidechainMC.next_resnum(sc_res_num)
	sidechainMC.apply(p)

	
#### BACKRUB PROTOCOL OPTIONS ####
prob_rotamer_only = 0.25
prob_rotamer_post_backrub = 0.25
prob_rotamer_pbr = 0.25

#### BACKRUB PROTOCOL ####
def Backrub_Mover(p):
	mc2.reset(p)
	for i in range(args.Number_of_Backrub_Moves):
		if rnd() < prob_rotamer_only:
			sidechainMC.apply(p)
		else:
			mc1.reset(p)
			backrub_mover.apply(p)
			mc1.boltzmann(p)
			if rnd() < prob_rotamer_post_backrub:
				backrub_sidechain_mover(p)
				if rnd() < prob_rotamer_pbr:
					backrub_sidechain_mover(p)
		mc2.boltzmann(p)			

### RESETTING TASKS AND MOVERS ###
def reset_pack_and_min(mutant_pose):
	task = standard_packer_task(mutant_pose)
	task.restrict_to_repacking()
	pack_mover.task(task)
	pyrosetta.rosetta.basic.options.set_integer_option('packing:ex1:level', 1)
	pyrosetta.rosetta.basic.options.set_integer_option('packing:ex2:level', 1)
	pyrosetta.rosetta.basic.options.set_integer_option('packing:ex3:level', 1)
	pack_mover.apply(p)
	minmover.apply(p)
	pyrosetta.rosetta.basic.options.set_integer_option('packing:ex1:level', 7)
	pyrosetta.rosetta.basic.options.set_integer_option('packing:ex1:level', 7)
	pyrosetta.rosetta.basic.options.set_integer_option('packing:ex1:level', 7)
	sidechainMC.set_task(task)
	mc1.reset(mutant_pose)
	mc2.reset(mutant_pose)
	# Import Global Variables
	global running_TRESPEFRET_SUM
	global running_TRESPEFRET_SUM_1
	global running_TRESPEFRET_SUM_2
	global running_EFRETKsq_SUM
	global running_EFRETKsq_SUM_1
	global running_EFRETKsq_SUM_2
	global running_EFRETKsq_SUM_W
	global running_EFRET_SUM
	global running_EFRET_SUM_1
	global running_EFRET_SUM_2
	global running_EFRET_SUM_W
	global running_TRESP_count
	# Update Global Variables
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

### FLUORESCENCE CALCULATORS ###
# Setting up Variables
QuantumYield = 0.0
JOverInt = 0.0
tau_d = 0.0 #seconds
## Variables for Counters
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

###  THE RUN PROTOCOL ###
for CaM_Mutant in simulation_list:
	p.assign(starting_pose)
	Mutate(p, CaM_Mutant[0], CaM_Mutant[1], CaM_Mutant[2], CaM_Mutant[3])
	reset_pack_and_min(p)
	mutant_start_pose.assign(p)
	for i in range(args.Number_of_Output_Structures):
		p.assign(mutant_start_pose)
		Backrub_Mover(p)
		running_TRESP_count += 1
		### For CNF and Tyr
		if CaM_Mutant[0] == 'CNF' or CaM_Mutant[0] == 'TYR':
			EFRET = compute_EFRET(p, CaM_Mutant[0], CaM_Mutant[4], CaM_Mutant[1], CaM_Mutant[3])
			running_EFRET_SUM += EFRET
			running_EFRET_AVG = running_EFRET_SUM/running_TRESP_count
			EFRETKsq, K_sq, R = compute_K2_EFRET(p, CaM_Mutant[0], CaM_Mutant[4], CaM_Mutant[1], CaM_Mutant[3])
			running_EFRETKsq_SUM += EFRETKsq
			running_EFRETKsq_AVG = running_EFRETKsq_SUM/running_TRESP_count
			TRESPEFRET = compute_TRESP_EFRET(p, CaM_Mutant[0], CaM_Mutant[4], CaM_Mutant[1], CaM_Mutant[3])
			running_TRESPEFRET_SUM += TRESPEFRET
			running_TRESP_AVG = running_TRESPEFRET_SUM/running_TRESP_count
			outf = open(str(CaM_Mutant[0]) + '_' + str(CaM_Mutant[1]) + '_' + str(CaM_Mutant[2]) + '_' + str(CaM_Mutant[3]) + '_' + str(args.Number_of_Input_Structure) + ".score", 'a')
			pdb_out = str(CaM_Mutant[0]) + '_' + str(CaM_Mutant[1]) + '_' + str(CaM_Mutant[2]) + '_' + str(CaM_Mutant[3]) + '_' + str(args.Number_of_Input_Structure) + '_' + str(i) + ".pdb"
			outf.write("%s\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (pdb_out, SF1(p), CA_rmsd(mutant_start_pose, p), R, K_sq, EFRET, EFRETKsq, TRESPEFRET, running_EFRET_AVG, running_EFRETKsq_AVG, running_TRESP_AVG))
			if args.Dump_PDBs and args.Dump_PDBs == True:
				p.dump_pdb(pdb_out)
			outf.close()
			if args.Show_Results and args.Show_Results == True:
				pmm.apply(p)
				print('Instantaneous Energy:' + str(SF1(p)))
				print('Average EFRET:' + str(running_EFRET_AVG))
				print('Average EFRETKsq:' + str(running_EFRETKsq_AVG))
				print('Average TRESP EFRET:' + str(running_TRESP_AVG))
		### For Trp		
		if CaM_Mutant[0] == 'TRP':
			EFRET_1, EFRET_2, EFRET_W = compute_EFRET(p, CaM_Mutant[0], CaM_Mutant[4], CaM_Mutant[1], CaM_Mutant[3])
			running_EFRET_SUM_1 += EFRET_1
			running_EFRET_AVG_1 = running_EFRET_SUM_1/running_TRESP_count
			running_EFRET_SUM_2 += EFRET_2
			running_EFRET_AVG_2 = running_EFRET_SUM_2/running_TRESP_count
			running_EFRET_SUM_W += EFRET_W
			running_EFRET_AVG_W = running_EFRET_SUM_W/running_TRESP_count
			EFRETKsq_1, EFRETKsq_2, EFRETKsq_W, K_sq_1, K_sq_2, K_sq_W, R_1, R_2, R_W = compute_K2_EFRET(p, CaM_Mutant[0], CaM_Mutant[4], CaM_Mutant[1], CaM_Mutant[3])
			running_EFRETKsq_SUM_1 += EFRETKsq_1
			running_EFRETKsq_AVG_1 = running_EFRETKsq_SUM_1/running_TRESP_count
			running_EFRETKsq_SUM_2 += EFRETKsq_2
			running_EFRETKsq_AVG_2 = running_EFRETKsq_SUM_2/running_TRESP_count
			running_EFRETKsq_SUM_W += EFRETKsq_W
			running_EFRETKsq_AVG_W = running_EFRETKsq_SUM_W/running_TRESP_count
			TRESPEFRET_1,TRESPEFRET_2 = compute_TRESP_EFRET(p, CaM_Mutant[0], CaM_Mutant[4], CaM_Mutant[1], CaM_Mutant[3])
			running_TRESPEFRET_SUM_1 += TRESPEFRET_1
			running_TRESPEFRET_SUM_2 += TRESPEFRET_2
			running_TRESP_AVG_1 = running_TRESPEFRET_SUM_1/running_TRESP_count
			running_TRESP_AVG_2 = running_TRESPEFRET_SUM_2/running_TRESP_count
			outf = open(str(CaM_Mutant[0]) + '_' + str(CaM_Mutant[1]) + '_' + str(CaM_Mutant[2]) + '_' + str(CaM_Mutant[3]) + '_' + str(args.Number_of_Input_Structure) + ".score", 'a')
			pdb_out = str(CaM_Mutant[0]) + '_' + str(CaM_Mutant[1]) + '_' + str(CaM_Mutant[2]) + '_' + str(CaM_Mutant[3]) + '_' + str(args.Number_of_Input_Structure) + '_' + str(i) + ".pdb"
			outf.write("%s\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (pdb_out, SF1(p), CA_rmsd(mutant_start_pose, p), R_1, R_2, R_W, K_sq_1, K_sq_2, K_sq_W, EFRET_1, EFRET_2, EFRET_W, EFRETKsq_1, EFRETKsq_2, EFRETKsq_W, TRESPEFRET_1, TRESPEFRET_2, running_EFRET_AVG_1, running_EFRET_AVG_2, running_EFRET_AVG_W, running_EFRETKsq_AVG_1, running_EFRETKsq_AVG_2, running_EFRETKsq_AVG_W, running_TRESP_AVG_1, running_TRESP_AVG_2))
			if args.Dump_PDBs and args.Dump_PDBs == True:
				p.dump_pdb(pdb_out)
			outf.close()
			if args.Show_Results and args.Show_Results == True:
				pmm.apply(p)
				print('Instantaneous Energy:' + str(SF1(p)))
				print('Average EFRET:' + str(running_EFRET_AVG_1) + ' ' + str(running_EFRET_AVG_2) + ' ' + str(running_EFRET_AVG_W))
				print('Average EFRETKsq:' + str(running_EFRETKsq_AVG_1) + ' ' + str(running_EFRETKsq_AVG_2) + ' ' + str(running_EFRETKsq_AVG_W))
				print('Average TRESP EFRET:' + str(running_TRESP_AVG_1) + ' ' + str(running_TRESP_AVG_2))
