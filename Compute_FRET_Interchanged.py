import pyrosetta
from math import *
from math import exp, log, pi, sqrt
from random import random as rnd
from random import randint
import numpy as np

# Definitions for Vector Processing
## Definition for Vector Averaging
def vectoravg(vec1,vec2):
	avgvec = pyrosetta.rosetta.numeric.xyzVector_double_t(((vec1.x+vec2.x)/2.0),((vec1.y+vec2.y)/2.0),((vec1.z+vec2.z)/2.0))
	return avgvec
## Definition for Weighted Vector Averaging
def weightedvectoravg(weight1,vec1,weight2,vec2):
	weightavgvec = pyrosetta.rosetta.numeric.xyzVector_double_t((((weight1*vec1.x)+(weight2*vec2.x))/(weight1+weight2)),(((weight1*vec1.y)+(weight2*vec2.y))/(weight1+weight2)),(((weight1*vec1.z)+(weight2*vec2.z))/(weight1+weight2)))
	return weightavgvec
	
# Definitions for Calculating EFRET
## Calculation of Average EFRET (Distance exclusive variable)
def compute_EFRET(pose_in, Donor_Fluor, Donor_QY, Donor_Res, Thio_Res):
	## Thioamide	
	ST = pose_in.residue(Thio_Res).xyz("ST")
	C = pose_in.residue(Thio_Res).xyz("C")
	N = pose_in.residue(Thio_Res+1).xyz("N")
	CNmp = vectoravg(C,N)
	## Donor Fluor Variables
	if Donor_Fluor	== 'CNF':
		QuantumYield = 0.11
		JOverInt = 7*10**12
		tau_d = 7.0*10**(-9)*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## Cyanophenylalanine (from Goldberg et. al)
		CE1 = pose_in.residue(Donor_Res).xyz("CE1")
		CE2 = pose_in.residue(Donor_Res).xyz("CE2")
		CEmp = vectoravg(CE1,CE2)
		R = (CEmp - CNmp).norm()
		EFRET = 1/(1+(R/R_0)**6)
		return EFRET
	if Donor_Fluor	== 'TYR':
		QuantumYield = 0.14
		JOverInt = 2.21*10**12
		tau_d = 3.4*10**(-9)*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## Tyrosine (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		CE1 = pose_in.residue(Donor_Res).xyz("CD1")
		CE2 = pose_in.residue(Donor_Res).xyz("CD2")
		CEmp = vectoravg(CE1,CE2)
		R = (CEmp - CNmp).norm()
		EFRET = 1/(1+(R/R_0)**6)
		return EFRET
	if Donor_Fluor	== 'TRP': 
		QuantumYield = 0.13
		JOverInt = 1.63*10**9
		tau_d = 3.1*10**(-9)*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## Tryptophan (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		### First Transition 4.8 M-1 x cm-1 x 10^(-3)
		NE = pose_in.residue(Donor_Res).xyz("NE1")
		CE3 = pose_in.residue(Donor_Res).xyz("CE3")
		CEmp_1 = vectoravg(NE,CE3)
		### Second Transition 0.8 M-1 x cm-1 x 10^(-3)
		CG = pose_in.residue(Donor_Res).xyz("CG")
		CZ2 = pose_in.residue(Donor_Res).xyz("CZ2")
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
	return K_sq
	
### Compute EFRET using KappaSquared from Dipoles
def compute_K2_EFRET(pose_in, Donor_Fluor, Donor_QY, Donor_Res, Thio_Res):
	## Thioamide	
	ST = pose_in.residue(Thio_Res).xyz("ST")
	C = pose_in.residue(Thio_Res).xyz("C")
	N = pose_in.residue(Thio_Res+1).xyz("N")
	CNmp = vectoravg(C,N)
	## Donor Fluor Variables
	if Donor_Fluor	== 'CNF': 
		QuantumYield = 0.11
		JOverInt = 7*10**12
		tau_d = 7.0*10**(-9)*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## Cyanophenylalanine (from Goldberg et. al)
		CE1 = pose_in.residue(Donor_Res).xyz("CE1")
		CE2 = pose_in.residue(Donor_Res).xyz("CE2")
		CEmp = vectoravg(CE1,CE2)
		R = (CEmp - CNmp).norm()
		K_sq_1 = compute_K2(ST,CNmp,CEmp,CE1)
		EFRETKsq_1 = 1/(1+(R/(R_0*((K_sq_1/(2/3))**(1/6))))**6)
		K_sq_2 = compute_K2(ST,CNmp,CEmp,CE2)
		EFRETKsq_2 = 1/(1+(R/(R_0*((K_sq_1/(2/3))**(1/6))))**6)
		K_sq = (K_sq_1 +K_sq_2) / 2
		EFRETKsq = (EFRETKsq_1 + EFRETKsq_2) / 2
		return EFRETKsq, K_sq, R
	if Donor_Fluor	== 'TYR': 
		QuantumYield = 0.14
		JOverInt = 2.21*10**12
		tau_d = 3.4*10**(-9)*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## Tyrosine (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		CD1 = pose_in.residue(Donor_Res).xyz("CD1")
		CD2 = pose_in.residue(Donor_Res).xyz("CD2")
		CDmp = vectoravg(CD1,CD2)
		R = (CDmp - CNmp).norm()
		K_sq = compute_K2(ST,CNmp,CD1,CD2)
		EFRETKsq = 1/(1+(R/(R_0*((K_sq/(2/3))**(1/6))))**6)
		return EFRETKsq, K_sq, R
	if Donor_Fluor	== 'TRP': 
		QuantumYield = 0.13
		JOverInt = 1.63*10**9
		tau_d = 3.1*10**(-9)*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## Tryptophan (from Antosiewicz et. al. Computation of the Dipole Moments of Proteins Biophysical Journal Vol. 69 1995 1344-1354)
		### First Transition 4.8 M-1 x cm-1 x 10^(-3)
		NE = pose_in.residue(Donor_Res).xyz("NE1")
		CE3 = pose_in.residue(Donor_Res).xyz("CE3")
		CEmp_1 = vectoravg(NE,CE3)
		### Second Transition 0.8 M-1 x cm-1 x 10^(-3)
		CG = pose_in.residue(Donor_Res).xyz("CG")
		CZ2 = pose_in.residue(Donor_Res).xyz("CZ2")
		CEmp_2 = vectoravg(CG,CZ2)
		### Combined Transition
		CEmp_W = weightedvectoravg(4.8,CEmp_1,0.8,CEmp_2)
		R_1 = (CEmp_1 - CNmp).norm()
		R_2 = (CEmp_2 - CNmp).norm()
		R_W = (CEmp_W - CNmp).norm()
		K_sq_1 = compute_K2(ST,CNmp,NE,CE3)
		K_sq_2 = compute_K2(ST,CNmp,CG,CZ2)
		K_sq_W = 4.8/(4.8+0.8)*K_sq_1 + 0.8/(4.8+0.8)*K_sq_2
		EFRETKsq_1 = 1/(1+(R_1/(R_0*((K_sq_1/(2/3))**(1/6))))**6)
		EFRETKsq_2 = 1/(1+(R_2/(R_0*((K_sq_2/(2/3))**(1/6))))**6)
		EFRETKsq_W = 1/(1+(R_W/(R_0*((K_sq_W/(2/3))**(1/6))))**6)
		return EFRETKsq_1, EFRETKsq_2, EFRETKsq_W, K_sq_1, K_sq_2, K_sq_W, R_1, R_2, R_W

## Calculation of TRESP EFRET (Excitonic Coupling)
def compute_TRESP_EFRET(pose_in, Donor_Fluor, Donor_QY, Donor_Res, Thio_Res):
	## Constants
	e_constant = 1.602*10**(-19) # Coulombs
	D_Cm_Conv = 3.3356*10**(-30) # Debye to Coulomb meter conversion
	## All listed charges are in Atomic Units from Multiwfn fits from Gaussian CIS-aug-cc-PVTZ simulations
	## Thioamide Transition Charges
	ST_xyz = np.array(pose_in.residue(Thio_Res).xyz("ST"))
	C_xyz = np.array(pose_in.residue(Thio_Res).xyz("C"))
	N_xyz = np.array(pose_in.residue(Thio_Res+1).xyz("N"))
	CA1_xyz = np.array(pose_in.residue(Thio_Res).xyz("CA"))
	CA2_xyz = np.array(pose_in.residue(Thio_Res+1).xyz("CA"))
	atom_i_list_xyz = [ST_xyz, C_xyz, N_xyz]#, CA1_xyz, CA2_xyz]	
	ST_chg = 0.001093/sqrt(2)#0.012159/sqrt(2)
	C_chg = -0.001682/sqrt(2)#-0.006775/sqrt(2)
	N_chg = 0.000589/sqrt(2)#-0.010651/sqrt(2)
	#CA1_chg = -0.007636/sqrt(2)
	#CA2_chg = 0.002535/sqrt(2)
	atom_i_list_chg = [ST_chg, C_chg, N_chg]#, CA1_chg, CA2_chg]
	if Donor_Fluor	== 'CNF':
		QuantumYield = 0.11
		JOverInt = 7*10**12
		tau_d = 7.0*10**(-9)#*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## CNF Transition Charges from S0S1 of Cyanobenzene
		CE1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE1"))
		CE2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE2"))
		CD1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD1"))
		CD2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD2"))
		CG_xyz = np.array(pose_in.residue(Donor_Res).xyz("CG"))
		CZ_xyz = np.array(pose_in.residue(Donor_Res).xyz("CZ"))
		CT_xyz = np.array(pose_in.residue(Donor_Res).xyz("CT"))
		NI_xyz = np.array(pose_in.residue(Donor_Res).xyz("NI"))
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
		TRESPEFRET_1 = 1/(1+(1/(tau_d*TRESP_k)))
		## COMPUTE A SECOND TIME WITH THE POSITIONS ACROSS THE SYSTEM FLIPPED TO SEE IF THERE IS AN IMPACT
				## CNF Transition Charges from S0S1 of Cyanobenzene
		CE2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE1"))
		CE1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE2"))
		CD2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD1"))
		CD1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD2"))
		CG_xyz = np.array(pose_in.residue(Donor_Res).xyz("CG"))
		CZ_xyz = np.array(pose_in.residue(Donor_Res).xyz("CZ"))
		CT_xyz = np.array(pose_in.residue(Donor_Res).xyz("CT"))
		NI_xyz = np.array(pose_in.residue(Donor_Res).xyz("NI"))
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
		TRESPEFRET_2 = 1/(1+(1/(tau_d*TRESP_k)))
		TRESPEFRET = (TRESPEFRET_1 + TRESPEFRET_2)/2
		return TRESPEFRET
	if Donor_Fluor	== 'TYR':
		QuantumYield = 0.14
		JOverInt = 2.21*10**12
		tau_d = 3.4*10**(-9)#*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## TYR Transition Charges from S0S1 of Methylphenol	
		CE2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE1"))
		CE1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE2"))
		CD2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD1"))
		CD1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD2"))
		CG_xyz = np.array(pose_in.residue(Donor_Res).xyz("CG"))
		CZ_xyz = np.array(pose_in.residue(Donor_Res).xyz("CZ"))
		OH_xyz = np.array(pose_in.residue(Donor_Res).xyz("OH"))
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
		TRESPEFRET_1 = 1/(1+(1/(tau_d*TRESP_k)))
		CE1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE1"))
		CE2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE2"))
		CD1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD1"))
		CD2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD2"))
		CG_xyz = np.array(pose_in.residue(Donor_Res).xyz("CG"))
		CZ_xyz = np.array(pose_in.residue(Donor_Res).xyz("CZ"))
		OH_xyz = np.array(pose_in.residue(Donor_Res).xyz("OH"))
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
		TRESPEFRET_2 = 1/(1+(1/(tau_d*TRESP_k)))
		TRESPEFRET = (TRESPEFRET_1 + TRESPEFRET_2)/2
		return TRESPEFRET
	if Donor_Fluor	== 'TRP':
		QuantumYield = 0.13
		JOverInt = 1.63*10**9
		tau_d = 3.1*10**(-9)#*Donor_QY/QuantumYield
		## Compute the R0
		R_0 = 0.211*((Donor_QY*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
		## TRP Transition Charges from S0S1 of Methyltryptophan
		CE2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE2"))
		CE3_xyz = np.array(pose_in.residue(Donor_Res).xyz("CE3"))
		CD1_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD1"))
		CD2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CD2"))
		CG_xyz = np.array(pose_in.residue(Donor_Res).xyz("CG"))
		CZ2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CZ2"))
		CZ3_xyz = np.array(pose_in.residue(Donor_Res).xyz("CZ3"))
		NE1_xyz = np.array(pose_in.residue(Donor_Res).xyz("NE1"))
		CB_xyz = np.array(pose_in.residue(Donor_Res).xyz("CB"))
		CH2_xyz = np.array(pose_in.residue(Donor_Res).xyz("CH2"))
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


