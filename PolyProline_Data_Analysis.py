import numpy as np
import sys
import glob
import argparse
from scipy import optimize

parser = argparse.ArgumentParser(description='Program')
parser.add_argument('-F', '--Donor_Fluor', action='store', type=str, required=True,
	help='3-Letter code for Donor fluorophore (CNF=cyanophenylalanine, TYR=tyrosine, TRP=tryptophan')
parser.add_argument('-D', '--Directory_Name', action='store', type=str, required=True,
	help='Name of the directory that houses outputs from all PolyProlines')
parser.add_argument('-S', '--Sampling_Name', action='store', type=str, required=False,
	help='Name of the sampling scheme used. If no scheme was specified do not use this flag')
args = parser.parse_args()

# Importing the correct files
import_set = []
if args.Sampling_Name:
	import_set = args.Directory_Name 
	if args.Directory_Name[len(args.Directory_Name)-1] == '/':
		import_set += '*/' + args.Donor_Fluor + '_*_' + args.Sampling_Name + '.score'
	else:	
		import_set += '*/' + args.Donor_Fluor + '_*_' + args.Sampling_Name + '.score'
	import_set = glob.glob(import_set)
	
else:
	if args.Directory_Name[len(args.Directory_Name)-1] == '/':
		import_set_1 = glob.glob(args.Directory_Name + '*/' + args.Donor_Fluor + '_*.score')
		import_set_2 = glob.glob(args.Directory_Name + '*/' + args.Donor_Fluor + '_*_*.score')
		import_set = list(set(import_set_1).difference(set(import_set_2)))
	else:	
		import_set_1 = glob.glob(args.Directory_Name + '/*/' + args.Donor_Fluor + '_*.score')
		import_set_2 = glob.glob(args.Directory_Name + '/*/' + args.Donor_Fluor + '_*_*.score')
		import_set = list(set(import_set_1).difference(set(import_set_2)))

# Setting Fluorophore Specifics
qeff = []
import_dtype = []
QuantumYield = 0
Lifetime = 0
JOverInt = 0

if args.Donor_Fluor == 'CNF':
	qeff = [0.698, 0.359, 0.221, 0.207, 0.066]
	import_dtype=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R', '<f8'), ('K_sq', '<f8'), ('ThetaD', '<f8'), ('ThetaA', '<f8'), ('ThetaDA', '<f8'), ('EFRET', '<f8'), ('EFRETKsq', '<f8'), ('TRESPEFRET', '<f8'), ('Run_EFRET_AVG', '<f8'), ('Run_EFRETKsq_AVG', '<f8'), ('Run_TRESPEFRET_AVG', '<f8')]
	QuantumYield = 0.11
	JOverInt = 7*10**12
	Lifetime = 7.0*10**(-9)
if args.Donor_Fluor == 'TYR':
	qeff = [0.914, 0.876, 0.705, 0.516, 0.362, 0.207, 0.107, -0.064, -0.090]
	import_dtype=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R', '<f8'), ('K_sq', '<f8'), ('ThetaD', '<f8'), ('ThetaA', '<f8'), ('ThetaDA', '<f8'), ('EFRET', '<f8'), ('EFRETKsq', '<f8'), ('TRESPEFRET', '<f8'), ('Run_EFRET_AVG', '<f8'), ('Run_EFRETKsq_AVG', '<f8'), ('Run_TRESPEFRET_AVG', '<f8')]
	QuantumYield = 0.14
	JOverInt = 2.21*10**12
	Lifetime = 3.4*10**(-9)
if args.Donor_Fluor == 'TRP':
	qeff = [0.688, 0.441, 0.150, 0.019, -0.038]	# 0.004
	import_dtype=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R_1', '<f8'), ('R_2', '<f8'), ('R_W', '<f8'), ('K_sq_1', '<f8'), ('ThetaD_1', '<f8'), ('ThetaA_1', '<f8'), ('ThetaDA_1', '<f8'), ('K_sq_2', '<f8'), ('ThetaD_2', '<f8'), ('ThetaA_2', '<f8'), ('ThetaDA_2', '<f8'), ('K_sq_W', '<f8'), ('EFRET_1', '<f8'), ('EFRET_2', '<f8'), ('EFRET_W', '<f8'), ('EFRETKsq_1', '<f8'), ('EFRETKsq_2', '<f8'), ('EFRETKsq_W', '<f8'), ('TRESPEFRET_1', '<f8'), ('TRESPEFRET_2', '<f8'), ('Run_EFRET_AVG_1', '<f8'), ('Run_EFRET_AVG_2', '<f8'), ('Run_EFRET_AVG_W', '<f8'), ('Run_EFRETKsq_AVG_1', '<f8'), ('Run_EFRETKsq_AVG_2', '<f8'), ('Run_EFRETKsq_AVG_W', '<f8'), ('Run_TRESPEFRET_AVG_1', '<f8'), ('Run_TRESPEFRET_AVG_2', '<f8')]
	JOverInt = 1.63*10**9
	QuantumYield = 0.13
	Lifetime = 3.1*10**(-9)

R_0 = 0.211*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
rate_R_0 = 9.78*10**(3)*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
rate_kf = QuantumYield/Lifetime
rate_knr = (1/Lifetime)-rate_kf

# Setting Up Alternative Quenching Models 
## Distance Dependent Quenching
ddqfit = lambda p, R_in: 1/(1+(p[0]*np.exp(((2*R_in) - p[1])/p[2])))
ddqerrfunc = lambda p, R_in, ef: ddqfit(p, R_in) - ef

## Distance Dependent Quenching with FRET
ddqfreterrfunc = lambda p, R_in, ef, eftresp: ddqfit(p, R_in) - ef + eftresp
#ddqfretfit = lambda p, R_in, Ksq_in: 1-(1/QuantumYield)*(rate_kf/(rate_kf + rate_knr + (p[0]*np.exp(-(R_in - p[1])/p[2])) + (1/Lifetime)*((R_0*((Ksq_in/(2/3))**(1/6)))/R_in)**6))
#ddqfreterrfunc = lambda p, R_in, Ksq_in, ef: ddqfretfit(p, R_in, Ksq_in) - ef
## Dexter Energy/Electron Transfer
dexfit = lambda p, R_in: 1/(1+(p[0]*np.exp(((2*R_in))/p[1])))
dexerrfunc = lambda p, R_in, ef: dexfit(p, R_in) - ef

## Dexter With FRET
dexfreterrfunc = lambda p, R_in, ef, eftresp: dexfit(p, R_in) - ef + eftresp

# Setting Up Data Holders
list_of_score_files = import_set
holder_dtype = []
ksq_holder_dtype = []
ksq_holder_dtype = [('P2_K_sq', '<f8'), ('P2_ThetaD', '<f8'), ('P2_ThetaA', '<f8'), ('P2_ThetaDA', '<f8'), ('P2_Distance', '<f8'), ('P3_K_sq', '<f8'), ('P3_ThetaD', '<f8'), ('P3_ThetaA', '<f8'), ('P3_ThetaDA', '<f8'), ('P3_Distance', '<f8'), ('P4_K_sq', '<f8'), ('P4_ThetaD', '<f8'), ('P4_ThetaA', '<f8'), ('P4_ThetaDA', '<f8'), ('P4_Distance', '<f8'), ('P5_K_sq', '<f8'), ('P5_ThetaD', '<f8'), ('P5_ThetaA', '<f8'), ('P5_ThetaDA', '<f8'), ('P5_Distance', '<f8'), ('P6_K_sq', '<f8'), ('P6_ThetaD', '<f8'), ('P6_ThetaA', '<f8'), ('P6_ThetaDA', '<f8'), ('P6_Distance', '<f8'), ('P7_K_sq', '<f8'), ('P7_ThetaD', '<f8'), ('P7_ThetaA', '<f8'), ('P7_ThetaDA', '<f8'), ('P7_Distance', '<f8'),  ('P8_K_sq', '<f8'), ('P8_ThetaD', '<f8'), ('P8_ThetaA', '<f8'), ('P8_ThetaDA', '<f8'), ('P8_Distance', '<f8'),  ('P9_K_sq', '<f8'), ('P9_ThetaD', '<f8'), ('P9_ThetaA', '<f8'), ('P9_ThetaDA', '<f8'),  ('P9_Distance', '<f8'),('P10_K_sq', '<f8'), ('P10_ThetaD', '<f8'), ('P10_ThetaA', '<f8'), ('P10_ThetaDA', '<f8'), ('P10_Distance', '<f8')]
ksqw_holder_dtype = [('P2_K_sq', '<f8'), ('P3_K_sq', '<f8'), ('P4_K_sq', '<f8'), ('P5_K_sq', '<f8'), ('P6_K_sq', '<f8'), ('P7_K_sq', '<f8'), ('P8_K_sq', '<f8'), ('P9_K_sq', '<f8'), ('P10_K_sq', '<f8')]

if args.Donor_Fluor == 'TRP':
	holder_dtype = [('Pro_Num', '<f8'), ('Exp_EFRET', '<f8'), ('Exp_EFRET_Std', '<f8'), ('Score', '<f8'), ('Score_Std', '<f8'), ('R', '<f8'), ('R_std', '<f8'), ('K_sq_1', '<f8'), ('K_sq_1_std', '<f8'), ('K_sq_2', '<f8'), ('K_sq_2_std', '<f8'), ('K_sq_W', '<f8'), ('K_sq_W_std', '<f8'), ('EFRET_1', '<f8'), ('EFRET_1_std', '<f8'), ('EFRET_2', '<f8'), ('EFRET_2_std', '<f8'), ('EFRET_W', '<f8'), ('EFRET_W_std', '<f8'), ('EFRETKsq_1', '<f8'), ('EFRETKsq_1_std', '<f8'), ('EFRETKsq_2', '<f8'), ('EFRETKsq_2_std', '<f8'), ('EFRETKsq_W', '<f8'), ('EFRETKsq_W_std', '<f8'), ('TRESPEFRET_1', '<f8'), ('TRESPEFRET_1_std', '<f8'), ('TRESPEFRET_2', '<f8'), ('TRESPEFRET_2_std', '<f8'), ('DDQ', '<f8'), ('DDQ_std', '<f8'), ('DEX', '<f8'), ('DEX_std', '<f8'), ('DDQ_FRET', '<f8'), ('DDQ_FRET_std', '<f8') , ('DDQFRET', '<f8'), ('DDQFRET_std', '<f8'), ('DEX_FRET', '<f8'), ('DEX_FRET_std', '<f8'), ('DEXFRET', '<f8'), ('DEXFRET_std', '<f8')]

else:
	holder_dtype = [('Pro_Num', '<f8'), ('Exp_EFRET', '<f8'), ('Exp_EFRET_Std', '<f8'), ('Score', '<f8'), ('Score_Std', '<f8'), ('R', '<f8'), ('R_std', '<f8'), ('K_sq', '<f8'), ('K_sq_std', '<f8'), ('EFRET', '<f8'), ('EFRET_std', '<f8'), ('EFRETKsq', '<f8'), ('EFRETKsq_std', '<f8'), ('TRESPEFRET', '<f8'), ('TRESPEFRET_std', '<f8'), ('DDQ', '<f8'), ('DDQ_std', '<f8'), ('DEX', '<f8'), ('DEX_std', '<f8'), ('DDQ_FRET', '<f8'), ('DDQ_FRET_std', '<f8') , ('DDQFRET', '<f8'), ('DDQFRET_std', '<f8'), ('DEX_FRET', '<f8'), ('DEX_FRET_std', '<f8'), ('DEXFRET', '<f8'), ('DEXFRET_std', '<f8')]

data_holder = np.zeros([5,], dtype=holder_dtype)
test_dataset = np.genfromtxt(list_of_score_files[0], dtype=import_dtype, delimiter='\t')
Ksq_holder = np.zeros([len(test_dataset),], dtype=ksq_holder_dtype)
Ksq_holder_2 = np.zeros([len(test_dataset),], dtype=ksq_holder_dtype)
Ksq_holder_w = np.zeros([len(test_dataset),], dtype=ksqw_holder_dtype)
	
# Performing Data Extraction
data_holder_set = []
for score_file in list_of_score_files:
	score_file_data = np.genfromtxt(score_file, dtype=import_dtype, delimiter='\t')
	proline_num = int(score_file_data[0][0].split('_')[1].split('P')[1])
	proline_idx = proline_num-2
	data_holder['Pro_Num'][proline_idx] = proline_num
	data_holder['Exp_EFRET'][proline_idx] = qeff[proline_idx]
	data_holder['Score'][proline_idx] = np.average(score_file_data['Score'])
	data_holder['Score_Std'][proline_idx] = np.std(score_file_data['Score'])
	if args.Donor_Fluor == 'TRP':
		data_holder['R'][proline_idx] = np.average(score_file_data['R_W'])
		data_holder['K_sq_1'][proline_idx] = np.average(score_file_data['K_sq_1'])
		data_holder['K_sq_2'][proline_idx] = np.average(score_file_data['K_sq_2'])
		data_holder['K_sq_W'][proline_idx] = np.average(score_file_data['K_sq_W'])
		data_holder['EFRET_1'][proline_idx] = np.average(score_file_data['EFRET_1'])
		data_holder['EFRET_2'][proline_idx] = np.average(score_file_data['EFRET_2'])
		data_holder['EFRET_W'][proline_idx] = np.average(score_file_data['EFRET_W'])
		data_holder['EFRETKsq_1'][proline_idx] = np.average(score_file_data['EFRETKsq_1'])
		data_holder['EFRETKsq_2'][proline_idx] = np.average(score_file_data['EFRETKsq_2'])
		data_holder['EFRETKsq_W'][proline_idx] = np.average(score_file_data['EFRETKsq_W'])
		data_holder['TRESPEFRET_1'][proline_idx] = np.average(score_file_data['TRESPEFRET_1'])
		data_holder['TRESPEFRET_2'][proline_idx] = np.average(score_file_data['TRESPEFRET_2'])
		data_holder['R_std'][proline_idx] = np.std(score_file_data['R_W'])
		data_holder['K_sq_1_std'][proline_idx] = np.std(score_file_data['K_sq_1'])
		data_holder['K_sq_2_std'][proline_idx] = np.std(score_file_data['K_sq_2'])
		data_holder['K_sq_W_std'][proline_idx] = np.std(score_file_data['K_sq_W'])
		data_holder['EFRET_1_std'][proline_idx] = np.std(score_file_data['EFRET_1'])
		data_holder['EFRET_2_std'][proline_idx] = np.std(score_file_data['EFRET_2'])
		data_holder['EFRET_W_std'][proline_idx] = np.std(score_file_data['EFRET_W'])
		data_holder['EFRETKsq_1_std'][proline_idx] = np.std(score_file_data['EFRETKsq_1'])
		data_holder['EFRETKsq_2_std'][proline_idx] = np.std(score_file_data['EFRETKsq_2'])
		data_holder['EFRETKsq_W_std'][proline_idx] = np.std(score_file_data['EFRETKsq_W'])
		data_holder['TRESPEFRET_1_std'][proline_idx] = np.std(score_file_data['TRESPEFRET_1'])	
		data_holder['TRESPEFRET_2_std'][proline_idx] = np.std(score_file_data['TRESPEFRET_2'])	
	else:	
		data_holder['R'][proline_idx] = np.average(score_file_data['R'])
		data_holder['K_sq'][proline_idx] = np.average(score_file_data['K_sq'])
		data_holder['EFRET'][proline_idx] = np.average(score_file_data['EFRET'])
		data_holder['EFRETKsq'][proline_idx] = np.average(score_file_data['EFRETKsq'])
		data_holder['TRESPEFRET'][proline_idx] = np.average(score_file_data['TRESPEFRET'])
		data_holder['R_std'][proline_idx] = np.std(score_file_data['R'])
		data_holder['K_sq_std'][proline_idx] = np.std(score_file_data['K_sq'])
		data_holder['EFRET_std'][proline_idx] = np.std(score_file_data['EFRET'])
		data_holder['EFRETKsq_std'][proline_idx] = np.std(score_file_data['EFRETKsq'])
		data_holder['TRESPEFRET_std'][proline_idx] = np.std(score_file_data['TRESPEFRET'])	
	## Pulling out the Kappa Squared Data
	if args.Donor_Fluor == 'TRP':
		Ksq_holder_w['P'+str(proline_idx+2)+'_K_sq'] = score_file_data['K_sq_W']
		Ksq_holder_2['P'+str(proline_idx+2)+'_K_sq'] = score_file_data['K_sq_2']
		Ksq_holder_2['P'+str(proline_idx+2)+'_ThetaD'] = score_file_data['ThetaD_2']
		Ksq_holder_2['P'+str(proline_idx+2)+'_ThetaA'] = score_file_data['ThetaA_2']
		Ksq_holder_2['P'+str(proline_idx+2)+'_ThetaDA'] = score_file_data['ThetaDA_2']
		Ksq_holder['P'+str(proline_idx+2)+'_K_sq'] = score_file_data['K_sq_1']
		Ksq_holder['P'+str(proline_idx+2)+'_ThetaD'] = score_file_data['ThetaD_1']
		Ksq_holder['P'+str(proline_idx+2)+'_ThetaA'] = score_file_data['ThetaA_1']
		Ksq_holder['P'+str(proline_idx+2)+'_ThetaDA'] = score_file_data['ThetaDA_1']
	else:
		Ksq_holder['P'+str(proline_idx+2)+'_K_sq'] = score_file_data['K_sq']
		Ksq_holder['P'+str(proline_idx+2)+'_ThetaD'] = score_file_data['ThetaD']
		Ksq_holder['P'+str(proline_idx+2)+'_ThetaA'] = score_file_data['ThetaA']
		Ksq_holder['P'+str(proline_idx+2)+'_ThetaDA'] = score_file_data['ThetaDA']
		Ksq_holder['P'+str(proline_idx+2)+'_Distance'] = score_file_data['R']
	
# Performing Alternative Transfer Mechanism Fitting
## Setting up initial conditions
pdex = [0.6, 7.1] # k, LDex
pddq = [0.6, 3.0, 5.0] # ka, a, re
pdexfret = [0.6, 7.1] # k, LDex
pddqfret = [0.6, 3.0, 5.0] # ka, a, re

input_R = data_holder['R']
if args.Donor_Fluor == 'TRP':
	#input_Ksq = data_holder['K_sq_1']
	input_eftresp = data_holder['TRESPEFRET_1']
else:
	#input_Ksq = data_holder['K_sq']
	input_eftresp = data_holder['TRESPEFRET']
input_qeff = data_holder['Exp_EFRET']


## Performing fitting
pddq1, success = optimize.leastsq(ddqerrfunc, pddq[:], args=(input_R, input_qeff), maxfev=16000)
data_holder['DDQ'] = ddqfit(pddq1, data_holder['R'])
pddqfret1, success = optimize.leastsq(ddqfreterrfunc, pddqfret[:], args=(input_R, input_qeff, input_eftresp), maxfev=16000)
data_holder['DDQ_FRET'] = ddqfit(pddqfret1, data_holder['R'])
data_holder['DDQFRET'] = data_holder['DDQ_FRET'] + input_eftresp
pdex1, success = optimize.leastsq(dexerrfunc, pdex[:], args=(input_R, input_qeff), maxfev=16000)
data_holder['DEX'] = dexfit(pdex1, data_holder['R'])
pdexfret1, success = optimize.leastsq(dexfreterrfunc, pdex[:], args=(input_R, input_qeff, input_eftresp), maxfev=16000)
data_holder['DEX_FRET'] = dexfit(pdexfret1, data_holder['R'])
data_holder['DEXFRET'] = data_holder['DEX_FRET'] + input_eftresp

# Formatting Outputs
## Output Naming
if args.Donor_Fluor == 'TRP':
	compiled_fret_data_header = 'Pro_Num \t Exp_EFRET \t Exp_EFRET_Std \t Score \t Score_Std \t R \t R_std \t K_sq_1 \t K_sq_1_std \t K_sq_2 \t K_sq_2_std \t K_sq_W \t K_sq_W_std \t EFRET_1 \t EFRET_1_std \t EFRET_2 \t EFRET_2_std \t EFRET_W \t EFRET_W_std \t EFRETKsq_1 \t EFRETKsq_1_std \t EFRETKsq_2 \t EFRETKsq_2_std \t EFRETKsq_W \t EFRETKsq_W_std \t TRESPEFRET_1 \t TRESPEFRET_1_std \t TRESPEFRET_2 \t TRESPEFRET_2_std \t DDQ \t DDQ_std \t DEX \t DEX_std \t DDQ_FRET \t DDQ_FRET_std \t DDQFRET \t DDQFRET_std \t DEX_FRET \t DEX_FRET_std \t DEXFRET \t DEXFRET_std'
else:	
	compiled_fret_data_header = 'Pro_Num \t Exp_EFRET \t Exp_EFRET_std \t Score \t Score_Std \t R \t R_std \t K_sq \t K_sq_std \t EFRET \t EFRET_std \t EFRETKsq \t EFRETKsq_std \t TRESPEFRET \t TRESPFRET_std \t DDQ \t DDQ_Std \t DEX \t DEX_Std \t DDQ_FRET \t DDQ_FRET_std \t DDQFRET \t DDQFRET_std \t DEX_FRET \t DEX_FRET_std \t DEXFRET \t DEXFRET_std'
kappa_squared_data_header = 'P2_K_sq \t P2_ThetaD \t P2_ThetaA \t P2_ThetaDA \t P2_Distance \t P3_K_sq \t P3_ThetaD \t P3_ThetaA \t P3_ThetaDA \t P3_Distance \t P4_K_sq \t P4_ThetaD \t P4_ThetaA \t P4_ThetaDA \t P4_Distance \t P5_K_sq \t P5_ThetaD \t P5_ThetaA \t P5_ThetaDA \t P5_Distance \t P6_K_sq \t P6_ThetaD \t P6_ThetaA \t P6_ThetaDA \t P6_Distance \t P7_K_sq \t P7_ThetaD \t P7_ThetaA \t P7_ThetaDA \t P7_Distance \t P8_K_sq \t P8_ThetaD \t P8_ThetaA \t P8_ThetaDA \t P8_Distance \t P9_K_sq \t P9_ThetaD \t P9_ThetaA \t P9_ThetaDA \t P9_Distance \t P10_K_sq \t P10_ThetaD \t P10_ThetaA \t P10_ThetaDA \t P10_Distance'
pdex_header = 'k \t LDex'
pddq_header = 'ka \t a \t re'
pddqfret_header = 'ka \t a \t re'
pdexfret_header = 'k \t LDex'

Ksq_list_filename = ''
Ksq_2_list_filename = ''
Ksq_W_list_filename = ''
compiled_fret_data_filename = ''                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
if args.Sampling_Name:
	Ksq_list_filename = 'Ksq_list_' + args.Donor_Fluor + '_' + args.Sampling_Name + '.txt'
	compiled_fret_data_filename = 'EFRET_Data_Compiled_' + args.Donor_Fluor + '_' + args.Sampling_Name + '.txt'
	ddq_filename = 'DDQ_TERMS_' + args.Donor_Fluor + '_' + args.Sampling_Name + '.txt'
	dex_filename = 'DEX_TERMS_' + args.Donor_Fluor + '_' + args.Sampling_Name + '.txt'
	ddqfret_filename = 'DDQFRET_TERMS_' + args.Donor_Fluor + '_' + args.Sampling_Name + '.txt'
	dexfret_filename = 'DEXFRET_TERMS_' + args.Donor_Fluor + '_' + args.Sampling_Name + '.txt'
else:	
	Ksq_list_filename = 'Ksq_list_' + args.Donor_Fluor + '.txt'
	compiled_fret_data_filename = 'EFRET_Data_Compiled_' + args.Donor_Fluor + '.txt'
	ddq_filename = 'DDQ_TERMS_' + args.Donor_Fluor + '.txt'
	dex_filename = 'DEX_TERMS_' + args.Donor_Fluor + '.txt'
	ddqfret_filename = 'DDQFRET_TERMS_' + args.Donor_Fluor + '.txt'
	dexfret_filename = 'DEXFRET_TERMS_' + args.Donor_Fluor + '.txt'

if args.Donor_Fluor == 'TRP':
	Ksq_list_filename = Ksq_list_filename[:len(Ksq_list_filename)-4]
	Ksq_2_list_filename = Ksq_list_filename + '_2.txt'
	Ksq_W_list_filename = Ksq_list_filename + '_W.txt'	
	Ksq_list_filename += '_1.txt'

	
## Output Writing
np.savetxt(Ksq_list_filename, Ksq_holder, fmt='%s', delimiter='\t', newline='\n', header=kappa_squared_data_header)
np.savetxt(compiled_fret_data_filename, data_holder, fmt='%s', delimiter='\t', newline='\n', header=compiled_fret_data_header)
np.savetxt(ddq_filename, pddq1, fmt='%s', delimiter='\t', newline='\n', header=pddq_header)
np.savetxt(dex_filename, pdex1, fmt='%s', delimiter='\t', newline='\n', header=pdex_header)
np.savetxt(ddqfret_filename, pddqfret1, fmt='%s', delimiter='\t', newline='\n', header=pddqfret_header)
np.savetxt(dexfret_filename, pdexfret1, fmt='%s', delimiter='\t', newline='\n', header=pdexfret_header)
if args.Donor_Fluor == 'TRP':
	np.savetxt(Ksq_2_list_filename, Ksq_holder_2, fmt='%s', delimiter='\t', newline='\n', header=kappa_squared_data_header)
	np.savetxt(Ksq_W_list_filename, Ksq_holder_w, fmt='%s', delimiter='\t', newline='\n', header=kappa_squared_data_header)
