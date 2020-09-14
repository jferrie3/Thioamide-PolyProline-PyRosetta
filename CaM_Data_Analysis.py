import numpy as np
import sys
import glob
import argparse
from scipy import optimize

parser = argparse.ArgumentParser(description='Program')
#parser.add_argument('-F', '--Donor_Fluor', action='store', type=str, required=True,
#	help='3-Letter code for Donor fluorophore (CNF=cyanophenylalanine, TYR=tyrosine, TRP=tryptophan')
parser.add_argument('-S', '--Score_File_Gen_Name', action='store', type=str, required=True,
	help='Folder Name')
args = parser.parse_args()

# Importing the correct files
import_set = args.Score_File_Gen_Name + '/*/' + '*.score'

# Setting Fluorophore Specifics
## Setting Up Placeholders
qeff = []
import_dtype = []
QuantumYield = 0
Lifetime = 0
JOverInt = 0
R_0 = 0
rate_R_0 = 0
rate_kf = 0
rate_knr = 0
pdex = [0.0, 0.0] # k, LDex
pddq = [0.0, 0.0, 0.0] # ka, a, re

def Fluor_Vals(Fluor_ID):
	QuantumYield = 0
	Lifetime = 0
	JOverInt = 0
	pdex = [0.0, 0.0] # k, LDex
	pddq = [0.0, 0.0, 0.0] # ka, a, re
	if Fluor_ID == 'CNF':
		import_dtype=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R', '<f8'), ('K_sq', '<f8'), ('EFRET', '<f8'), ('EFRETKsq', '<f8'), ('TRESPEFRET', '<f8'), ('Run_EFRET_AVG', '<f8'), ('Run_EFRETKsq_AVG', '<f8'), ('Run_TRESPEFRET_AVG', '<f8')]
		QuantumYield = 0.11
		JOverInt = 7*10**12
		Lifetime = 7.0*10**(-9)
		pdex = [0.009329236720641823, 6.673761976674292] # k, LDex
		pddq = [0.12872527819891197, 17.515471468774475, 6.673763651770554] # ka, a, re
	if Fluor_ID == 'TYR':
		import_dtype=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R', '<f8'), ('K_sq', '<f8'), ('EFRET', '<f8'), ('EFRETKsq', '<f8'), ('TRESPEFRET', '<f8'), ('Run_EFRET_AVG', '<f8'), ('Run_EFRETKsq_AVG', '<f8'), ('Run_TRESPEFRET_AVG', '<f8')]
		QuantumYield = 0.14
		JOverInt = 2.21*10**12
		Lifetime = 3.4*10**(-9)
		pdex = [0.010563577334917673, 6.340720464017725] # k, LDex
		pddq = [17.672581972421234, 47.063076078015904, 6.3407117209127275] # ka, a, re
	if Fluor_ID == 'TRP':
		import_dtype=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R_1', '<f8'), ('R_2', '<f8'), ('R_W', '<f8'), ('K_sq_1', '<f8'), ('K_sq_2', '<f8'), ('K_sq_W', '<f8'), ('EFRET_1', '<f8'), ('EFRET_2', '<f8'), ('EFRET_W', '<f8'), ('EFRETKsq_1', '<f8'), ('EFRETKsq_2', '<f8'), ('EFRETKsq_W', '<f8'), ('TRESPEFRET_1', '<f8'), ('TRESPEFRET_2', '<f8'), ('Run_EFRET_AVG_1', '<f8'), ('Run_EFRET_AVG_2', '<f8'), ('Run_EFRET_AVG_W', '<f8'), ('Run_EFRETKsq_AVG_1', '<f8'), ('Run_EFRETKsq_AVG_2', '<f8'), ('Run_EFRETKsq_AVG_W', '<f8'), ('Run_TRESPEFRET_AVG_1', '<f8'), ('Run_TRESPEFRET_AVG_2', '<f8')]
		JOverInt = 1.63*10**9
		QuantumYield = 0.13
		Lifetime = 3.1*10**(-9)
		pdex = [0.011547369870029436, 6.136179129282161] # k, LDex
		pddq = [0.05978487794956127, 10.089823332063219, 6.136167071383317] # ka, a, re
	R_0 = 0.211*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
	rate_R_0 = 9.78*10**(3)*((QuantumYield*(2/3)*(1.33**(-4))*JOverInt))**(1/6)
	rate_kf = QuantumYield/Lifetime
	rate_knr = (1/Lifetime)-rate_kf
	return import_dtype, Lifetime, JOverInt, R_0, rate_R_0, rate_kf, rate_knr, pddq, pdex, QuantumYield

# Managing Different Labeling Locations
# Storage Format: Idx, QY, EFRET
label_site_dict = {'CNF_12_TBF_153':[0, 0.003, 0.11],
	'CNF_12_TBL_163':[1, 0.003, 0.37],
	'CNF_92_TBF_153':[2, 0.014, 0.79],
	'CNF_92_TBL_163':[3, 0.014, 0.43],
	'CNF_99_TBF_153':[4, 0.011, 0.25],
	'CNF_99_TBL_163':[5, 0.011, 0.17],
	'CNF_99_TBF_168':[6, 0.011, 0.10],
	'TYR_138_TBF_153':[7, 0.061, 0.50],
	'TYR_99_TBF_153':[9, 0.105, 0.22],
	'TRP_99_TBF_153':[8, 0.15, 0.27],
'TRP_138_TBF_153':[10, 0.29, 0.58]}

# Storing Imported Data
import_dtype_CNF_TYR=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R', '<f8'), ('K_sq', '<f8'), ('EFRET', '<f8'), ('EFRETKsq', '<f8'), ('TRESPEFRET', '<f8'), ('Run_EFRET_AVG', '<f8'), ('Run_EFRETKsq_AVG', '<f8'), ('Run_TRESPEFRET_AVG', '<f8')]
import_dtype_TRP=[('PDB_OUT', '<U18'), ('Score', '<f8'), ('CA_RMSD', '<f8'), ('R_1', '<f8'), ('R_2', '<f8'), ('R_W', '<f8'), ('K_sq_1', '<f8'), ('K_sq_2', '<f8'), ('K_sq_W', '<f8'), ('EFRET_1', '<f8'), ('EFRET_2', '<f8'), ('EFRET_W', '<f8'), ('EFRETKsq_1', '<f8'), ('EFRETKsq_2', '<f8'), ('EFRETKsq_W', '<f8'), ('TRESPEFRET_1', '<f8'), ('TRESPEFRET_2', '<f8'), ('Run_EFRET_AVG_1', '<f8'), ('Run_EFRET_AVG_2', '<f8'), ('Run_EFRET_AVG_W', '<f8'), ('Run_EFRETKsq_AVG_1', '<f8'), ('Run_EFRETKsq_AVG_2', '<f8'), ('Run_EFRETKsq_AVG_W', '<f8'), ('Run_TRESPEFRET_AVG_1', '<f8'), ('Run_TRESPEFRET_AVG_2', '<f8')]
import_data_CNF_12_TBF_153 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_CNF_12_TBL_163 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_CNF_92_TBF_153 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_CNF_92_TBL_163 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_CNF_99_TBF_153 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_CNF_99_TBL_163 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_CNF_99_TBF_168 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_TYR_138_TBF_153 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_TRP_99_TBF_153 = np.zeros([0,], dtype=import_dtype_TRP)
import_data_TYR_99_TBF_153 = np.zeros([0,], dtype=import_dtype_CNF_TYR)
import_data_TRP_138_TBF_153 = np.zeros([0,], dtype=import_dtype_TRP)
import_data_set = [import_data_CNF_12_TBF_153, import_data_CNF_12_TBL_163, import_data_CNF_92_TBF_153, import_data_CNF_92_TBL_163, import_data_CNF_99_TBF_153, import_data_CNF_99_TBL_163, import_data_CNF_99_TBF_168, import_data_TYR_138_TBF_153, import_data_TRP_99_TBF_153, import_data_TYR_99_TBF_153, import_data_TRP_138_TBF_153]

# Setting Up Alternative Quenching Models 
## Distance Dependent Quenching
ddqfit = lambda p, R_in, QY1, QY2: 1/(1+(QY1/QY2)*(p[0]*np.exp(((2*R_in) - p[1])/p[2])))
## Distance Dependent Quenching with FRET
## Dexter Energy/Electron Transfer
dexfit = lambda p, R_in: 1/(1+(p[0]*np.exp(((2*R_in))/p[1])))

# Setting Up Data Holders
holder_dtype = [('Label_Sites', '<U18'), ('XHOLDER', '<f8'), ('Exp_EFRET', '<f8'), ('Exp_EFRET_STD', '<f8'), ('Score', '<f8'), ('Score_STD', '<f8'), ('R', '<f8'), ('R_STD', '<f8'), ('K_sq', '<f8'), ('K_sq_STD', '<f8'), ('EFRET', '<f8'), ('EFRET_STD', '<f8'), ('EFRETKsq', '<f8'), ('EFRETKsq_STD', '<f8'), ('TRESPEFRET', '<f8'), ('TRESPEFRET_STD', '<f8'), ('DDQ', '<f8'), ('DDQ_STD', '<f8'), ('DEX', '<f8'), ('DEX_STD', '<f8')]
data_holder = np.zeros([11,], dtype=holder_dtype)
# Performing Data Extraction
data_holder_set = []
list_of_score_files = glob.glob(import_set)
for score_file in list_of_score_files:
	core_file_name_set = score_file.split('/')[len(score_file.split('/'))-1].split('_')
	core_file_name = core_file_name_set[0] + '_' + core_file_name_set[1] + '_' + core_file_name_set[2] + '_' + core_file_name_set[3]
	core_file_info_set = label_site_dict[core_file_name]
	import_dtype, Lifetime, JOverInt, R_0, rate_R_0, rate_kf, rate_knr, pddq, pdex, QuantumYield = Fluor_Vals(core_file_name_set[0]) 
	score_file_data = np.genfromtxt(score_file, dtype=import_dtype, delimiter='\t')
	import_data_set[label_site_dict[core_file_name][0]] = np.append(import_data_set[label_site_dict[core_file_name][0]], score_file_data)

for label_idx,import_data_item in enumerate(import_data_set):
	score_file = import_data_item['PDB_OUT'][0]
	core_file_name_set = score_file.split('/')[len(score_file.split('/'))-1].split('_')
	core_file_name = core_file_name_set[0] + '_' + core_file_name_set[1] + '_' + core_file_name_set[2] + '_' + core_file_name_set[3]
	data_holder['Label_Sites'][label_idx] = core_file_name
	data_holder['Exp_EFRET'][label_idx] = label_site_dict[core_file_name][2]
	data_holder['Score'][label_idx] = np.average(import_data_item['Score'])
	import_dtype, Lifetime, JOverInt, R_0, rate_R_0, rate_kf, rate_knr, pddq, pdex, QuantumYield = Fluor_Vals(core_file_name_set[0]) 
	if core_file_name_set[0] == 'CNF' or core_file_name_set[0] == 'TYR':
		import_data_item['Run_EFRET_AVG'] = dexfit(pdex, import_data_item['R'])
		import_data_item['Run_EFRETKsq_AVG'] = ddqfit(pddq, import_data_item['R'], QuantumYield, label_site_dict[core_file_name][1])
		data_holder['R'][label_idx] = np.average(import_data_item['R'])
		data_holder['K_sq'][label_idx] = np.average(import_data_item['K_sq'])
		data_holder['EFRET'][label_idx] = np.average(import_data_item['EFRET'])
		data_holder['EFRETKsq'][label_idx] = np.average(import_data_item['EFRETKsq'])
		data_holder['TRESPEFRET'][label_idx] = np.average(import_data_item['TRESPEFRET'])
		data_holder['DEX'][label_idx] = np.average(import_data_item['Run_EFRET_AVG'])
		data_holder['DDQ'][label_idx] = np.average(import_data_item['Run_EFRETKsq_AVG'])
		data_holder['R_STD'][label_idx] = np.std(import_data_item['R'])
		data_holder['K_sq_STD'][label_idx] = np.std(import_data_item['K_sq'])
		data_holder['EFRET_STD'][label_idx] = np.std(import_data_item['EFRET'])
		data_holder['EFRETKsq_STD'][label_idx] = np.std(import_data_item['EFRETKsq'])
		data_holder['TRESPEFRET_STD'][label_idx] = np.std(import_data_item['TRESPEFRET'])
		data_holder['DEX_STD'][label_idx] = np.std(import_data_item['Run_EFRET_AVG'])
		data_holder['DDQ_STD'][label_idx] = np.std(import_data_item['Run_EFRETKsq_AVG'])
	if core_file_name_set[0] == 'TRP':
		import_data_item['Run_EFRET_AVG_1'] = dexfit(pdex, import_data_item['R_1'])
		import_data_item['Run_EFRETKsq_AVG_1'] = ddqfit(pddq, import_data_item['R_1'], QuantumYield, label_site_dict[core_file_name][1])
		data_holder['R'][label_idx] = np.average(import_data_item['R_1'])
		data_holder['K_sq'][label_idx] = np.average(import_data_item['K_sq_1'])
		data_holder['EFRET'][label_idx] = np.average(import_data_item['EFRET_1'])
		data_holder['EFRETKsq'][label_idx] = np.average(import_data_item['EFRETKsq_1'])
		data_holder['TRESPEFRET'][label_idx] = np.average(import_data_item['TRESPEFRET_1'])
		data_holder['DEX'][label_idx] = np.average(import_data_item['Run_EFRET_AVG_1'])
		data_holder['DDQ'][label_idx] = np.average(import_data_item['Run_EFRETKsq_AVG_1'])
		data_holder['R_STD'][label_idx] = np.std(import_data_item['R_1'])
		data_holder['K_sq_STD'][label_idx] = np.std(import_data_item['K_sq_1'])
		data_holder['EFRET_STD'][label_idx] = np.std(import_data_item['EFRET_1'])
		data_holder['EFRETKsq_STD'][label_idx] = np.std(import_data_item['EFRETKsq_1'])
		data_holder['TRESPEFRET_STD'][label_idx] = np.std(import_data_item['TRESPEFRET_1'])
		data_holder['DEX_STD'][label_idx] = np.std(import_data_item['Run_EFRET_AVG_1'])
		data_holder['DDQ_STD'][label_idx] = np.std(import_data_item['Run_EFRETKsq_AVG_1'])

data_holder_header = 'Label_Sites \t XHOLDER \t Exp_EFRET \t Exp_EFRET_STD \t Score \t Score_STD \t R \t R_STD \t K_sq \t K_sq_STD \t EFRET \t EFRET_STD \t EFRETKsq \t EFRETKsq_STD \t TRESPEFRET \t TRESPEFRET_STD \t DDQ \t DDQ_STD \t DEX'
compiled_fret_data_filename = 'Analysis_of_' + args.Score_File_Gen_Name.split('/')[0] + '.txt'
np.savetxt(compiled_fret_data_filename, data_holder, fmt='%s', delimiter='\t', newline='\n', header=data_holder_header)

'''# Performing Alternative Transfer Mechanism Fitting
## Setting up initial conditions
pdex = [0.6, 7.1] # k, LDex
pddq = [0.6, 3.0, 5.0] # ka, a, re
pddqfret = [0.6, 3.0, 5.0] # ka, a, re
input_R = data_holder['R']
input_Ksq = data_holder['K_sq']
input_qeff = data_holder['Exp_EFRET']

## Performing fitting
pddq1, success = optimize.leastsq(ddqerrfunc, pddq[:], args=(input_R, input_qeff), maxfev=16000)
data_holder['DDQ'] = ddqfit(pddq1, data_holder['R'])
pddqfret1, success = optimize.leastsq(ddqfreterrfunc, pddqfret[:], args=(input_R, input_Ksq, input_qeff), maxfev=16000)
data_holder['DDQFRET'] = ddqfretfit(pddqfret, data_holder['R'], data_holder['K_sq'])
pdex1, success = optimize.leastsq(dexerrfunc, pdex[:], args=(input_R, input_qeff), maxfev=16000)
data_holder['DEX'] = dexfit(pddq1, data_holder['R'])
EFRET_RMSD = np.sqrt(np.average(np.square(data_holder['Exp_EFRET']-data_holder['EFRET'])))
'''