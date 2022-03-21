""" --- Notes ---

This script is for analysis of Allosteric signalling map (ASM) raw data obtained from AlloSigMA.
It was developed using Python 3.9.

    --Cautions--
     This default setting of this script is to perform analysis at per-residue resolution (one residue at a time)
     For anyone who wants to analyse effects of perturbating two or more residues at the same time,
     please define residues in mt_r_list manually:
     for example, to analyse effects of perturbating 230 and 231 (experiment 1), and 250 and 251 (experiment 2):
     mt_r_list = [['230', '231'], ['250', '251']]

    --Data--
     Materials:
     ASM raw data

     Require inputs:
     1. paths for data/outputs;
     2. Specific residues to look at (e.g. residue 383, residue 466);
     3. Specific regions to look at (e.g. Epoxide positioners (383 and 466))

     Output:
     1.json_all_singular.txt (Effects (delta H) of perturbating one single residue on specific residues;
     2.json_all.txt (Avaerage effects (delta H) of perturbating one signle residue on specific regions.

"""

import os
import json
import operator
import sys


#### --- Inputs --- ####
# Set Paths
current_path = sys.path[0]
os.chdir(current_path)

# Inputs 1: paths for data/outputs
ASM_data_path = '../raw_data/AlloSigMA/ASM/87TT043W/WorkFiles/JSONObjs'
output_path = '../outputs/ASM/'

ml_path = os.path.join(ASM_data_path,'apoCTDsEH_chainA_mutationlist.json')
rl_path = os.path.join(ASM_data_path,'apoCTDsEH_chainA_responselist.json')
up_path = os.path.join(ASM_data_path,'apoCTDsEH_chainA_ASM_UP.json')
dw_path = os.path.join(ASM_data_path,'apoCTDsEH_chainA_ASM_DOWN.json')

# Define protein residue id
start_resi = 230
end_resi = 546

# Inputs 2: Specific residues to look at
see_resi_list = ['383', '466',
                 '335', '496', '524',
                 '339', '499', '336', '469',
                 '417', '525', '267', '408', '419', '498']

# Inputs 2: Specific regions to look at
regions_list = [['383', '466'],
                ['335', '496', '524'],
                ['339', '499', '336', '469'],
                ['417', '525', '267', '408', '419', '498']]
regions_name = ['Epoxide positioners',
                'Catalytic Triad',
                'W336 niche',
                'F267 pocket']


# Generate list of perturbating residues for each experiment (default: running all residues one-by-one)
mt_r_list = []
for i in range(start_resi, end_resi):
    r_list = []
    r_list.append(str(i))
    mt_r_list.append(r_list)


#### --- Main --- ####
### Read raw_data ###
with open(ml_path, 'r') as jf:
    muta_l = json.load(jf)
    # print(muta_l)

with open(rl_path, 'r') as jf:
    resp_l = json.load(jf)
    # print(resp_l)

with open(up_path, 'r') as jf:
    up_l = json.load(jf)
    # print(up_l[0])

with open(dw_path, 'r') as jf:
    dw_l = json.load(jf)
    # print(dw_l[0])

### Data processing ###
## Calculate deltaH modulation ranges ##
# Step 1: calculate moduration range (delta Gi)
overall = []  # overall records modulation range (deltaGi)

for u, d in zip(up_l, dw_l):
    one_resi = []
    for ur, dr in zip(u, d):
        modu_range = ur - dr
        one_resi.append(modu_range)
    overall.append(one_resi)

# Step 2: calcultae moduration range (delta H)
modu_oa = []  # modu_oa records modulation range (delta H)
for il in overall:
    modu_oner = []
    mean = sum(il) / len(il)
    for i in il:
        mi = i - mean
        modu_oner.append(mi)
    modu_oa.append(modu_oner)

### Tidy up informations  ###
## for effects on specific residues entered as inputs 2
muta_all = []
response_all = []
value_all = []

for mt_r in mt_r_list:
    for see_resi in see_resi_list:
        muta_current = []
        # specify one resi wanted to check
        # specify looking sites and mutated residues
        posi_see = int(see_resi) - start_resi

        # reset b_factor to 0
        bf = []
        for i in range(0, len(resp_l)):
            bf.append(float(0))

        # calculate new b_factor
        all_posi = []
        for r in mt_r:
            mt_r_posi = [i for i, x in enumerate(muta_l) if (r in str(x))][0]
            all_posi.append(mt_r_posi)
            muta_current.append(r)

        for p in all_posi:
            add_l = modu_oa[p]
            bf = list(map(operator.add, bf, add_l))

        muta_all.append(muta_current)
        response_all.append(see_resi)
        value_all.append(bf[posi_see])

## for effects on specific regions entered as inputs 3
rg_fn_s = []
rg_fn_v = []
rg_fn_n = []
for regions, name in zip(regions_list, regions_name):
    rg_m = []
    rg_r = []
    rg_v = []

    for m, r, v in zip(muta_all, response_all, value_all):
        for rg in regions:
            if rg is r:
                rg_m.append(m[0])
                rg_r.append(r)
                rg_v.append(v)


    for mi in range(0, len(set(rg_m))):
        cur_position = [i for i, x in enumerate(rg_m) if x == list(set(rg_m))[mi]]
        rg_fn_s.append(list(set(rg_m))[mi])
        fn_v = []
        for cp in cur_position:
            fn_v.append(rg_v[cp])
        rg_fn_v.append(sum(fn_v) / len(fn_v))
        rg_fn_n.append(name)



#### --- Outputs --- ####
with open(os.path.join(output_path,'json_all_singular.txt'), 'w') as wf:
    wf.write('mutation;response;value \n')
    for m, re, v in zip(muta_all, response_all, value_all):
        wf.write("{0};{1};{2} \n".format(m, re, v))
    wf.close()

with open(os.path.join(output_path,'json_all.txt'), 'w') as wf:
    wf.write('mutation;position;value \n')
    for m, n, v in zip(rg_fn_s, rg_fn_n, rg_fn_v):
        wf.write("{0};{1};{2} \n".format(m, n, v))
    wf.close()

print('Finish!')
#END