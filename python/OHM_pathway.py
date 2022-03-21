""" --- Notes ---

This script is for analysis and visualisation of OHM pathway data.
It was developed using Python 3.9.
********* Please run this script in PyMOL. ********

    --Cautions--
    1. This script provides two types of analysis/visualisation methods:
        # state = 0: display the longest path of top5 pathways;
        # state = 1: display the top1 pathway.
    2. Codes of customize pathway cylinder displays are available at the end of this script.

    --Data--
     Materials:
     OHM data
     center_of_mass.py

     Require inputs:
     1. Specify path to OHM data.
     2. Specify state of analysis/visualisation.
     3. If customizing pathway cylinder displays, specify all pathway residues.

     Output:
     PyMOL windows with selected display.

"""



import os
import re
import sys
import center_of_mass
import pymol
from pymol.cgo import *

#### --- Inputs --- ####

# Inputs: path to OHM data
OHM_data_path = '../raw_data/OHM/506pathway'
os.chdir(OHM_data_path)
# Define state of display
# state = 0: display the longest path of top5 pathways;
# state = 1: display the top1 pathway

state = 1

one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
              'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
              'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
              'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

#### --- Main --- ####
all_files = os.listdir()

for f in all_files:
    if f.endswith('.pse'):
        ohm_pse_path = os.path.join(os.getcwd(), f)
    if f.endswith('.txt'):
        ohm_txt_path = os.path.join(os.getcwd(), f)


# 1. Obtain information from ranked pathway txt file
txt_info_ori = [] # all information
txt_otre = [] # information of pathway residues

with open(ohm_txt_path, 'r') as f:
    lines = f.readlines()
    for l in lines:
        l = l.replace('\ufeff', '').replace('\n', '')
        l_info = l.split(' ', maxsplit=2)
        l_otre = re.findall(r'\d+', l_info[-1])
        txt_info_ori.append(l_info)
        txt_otre.append(l_otre)


# 2. Open pse and obtain information
pymol.finish_launching()
pymol.cmd.load(ohm_pse_path)

# detect all objects+sele+path
name_lists = pymol.cmd.get_names('all')

# detect inputted residues
ip_f = [] # inputted functional residues (end points of pathway)
ip_a = [] # inputted allosteric residue (start points of pathway)

pymol.cmd.iterate('active', 'ip_f.append(resi)')
pymol.cmd.iterate('allo', 'ip_a.append(resi)')
ip_f = list(set(ip_f))
ip_a = list(set(ip_a))
print('Functional sites:', ip_f)
print('Allosteric sites:', ip_a)
# obtain residues in each pathway
py_p_all = [] # all residues
py_p_ot = [] # except inputted residues
py_p_ip = [] # inputted residues
for name in name_lists:
    if str(name).startswith('path') & (('residues') not in str(name)):  # equal to obtain path info
        n_path_sele = name + '-residues'
        p_all = []
        p_ot = []
        p_ip = []
        pymol.cmd.iterate(n_path_sele, 'p_all.append(resi)')
        p_all = list(set(p_all))
        py_p_all.append(p_all)
        for pa in p_all:
            if (pa not in ip_f) and (pa not in ip_a):
                p_ot.append(pa)
            else:
                p_ip.append(pa)
        py_p_ot.append(p_ot)
        py_p_ip.append(p_ip)
    if str(name).endswith('bfactor'):
        pro_name = name


# allocate pathways number in pse file to ranks, by:
# comparing py_p_ot (pathway residues in pse file) and txt_otre (pathway residues in txt files)
txt_p_match_to = []
for tot in txt_otre:
    mt = []
    mt = [ni for ni, x in enumerate(py_p_ot) if (sorted(x) == sorted(tot))]
    txt_p_match_to.append(mt)

print('matching information:')
print(txt_p_match_to)

# 3. display selected pathways (according to state parameters) in pymol
if state == 1:
    enable_ind = txt_p_match_to[0]
    enable_path = []
    for i in enable_ind:
        ei = int(i) + 1
        ep = 'path' + str(ei)
        enable_path.append(ep)
    print('*** Displaying pathway rank 0 in txt file *** ')
    print("Enable pse. pathways:")
    print(enable_path)

elif state == 0:  # select the longest path
    # TODO look for top5 paths with matching
    matching_pos = [i for i, x in enumerate(txt_p_match_to) if (bool(x))]
    top5_mp = matching_pos[0:5]
    top5_length = []
    for mp in top5_mp:
        top5_length.append(len(txt_otre[mp]))

    longest_pos = [a for a, b in zip(top5_mp, top5_length) if (b == max(top5_length))]
    if (len(longest_pos) > 1):
        fn_longest_pos = min(longest_pos)
    else:
        fn_longest_pos = longest_pos[0]
    # extract the longest path
    enable_ind = txt_p_match_to[fn_longest_pos]
    enable_path = []
    for i in enable_ind:
        ei = int(i) + 1
        ep = 'path' + str(ei)
        enable_path.append(ep)
    print('*** Displaying pathway rank ' + str(fn_longest_pos) + ' in txt file ***')
    print("Enable pse. pathways:")
    print(enable_path)
else:
    print('*********Cannot identify state**********')

# 4. Tidy up protein image, add center spheres to represent residues

pymol.viewing.label('all','')
pymol.viewing.disable('all')
pymol.viewing.enable(pro_name)
pymol.viewing.color('white', pro_name)
pymol.viewing.enable('active')
pymol.viewing.enable('allos')



for ep in enable_path:
    enable_residues = ep + '-residues and n.CA'
    eb_r = ep + '-residues'
    pymol.viewing.enable(ep)
    pymol.cmd.select(eb_r)
    pymol.viewing.label("sele and n. CA", '"%s%s" %(one_letter[resn],resi)')

# select all ot residues, generate pseudoatom and color with b-factors
if state == 1:
    ot_resi_lists = txt_otre[0]
elif state == 0:
    ot_resi_lists = txt_otre[fn_longest_pos]

for r in ot_resi_lists:
    # print(r)
    rn = '/' + pro_name + '///' + r
    pymol.cmd.select(rn)
    pymol.viewing.label("sele and n. CA", '"%s%s" %(one_letter[resn],resi)')
    center_of_mass.com('sele', object=r)
    pymol.viewing.color('yellow', r)

for f in ip_f:
    rn = '/active///' + f
    pymol.cmd.select(rn)
    # pymol.viewing.label('sele and n. CA', 'resi')
    center_of_mass.com('sele', object=f)
    pymol.viewing.color('marine', f)

for a in ip_a:
    rn = '/allos///' + a
    pymol.cmd.select(rn)
    # pymol.viewing.label('sele and n. CA', 'resi')
    center_of_mass.com('sele', object=a)
    pymol.editing.alter(a, "vdw=1.4")
    pymol.viewing.show('spheres', a)
    pymol.viewing.color('red', a)

# color functional sites
ep_site = '/active///383+466'
ct_site = '/active///335+496+524'
f267_site = '/active///417+525+267+408+419'
w336_site = '/active///339+499+336+469'

# pymol.cmd.select(vertex_site)
pymol.viewing.show('mesh', ep_site)
pymol.viewing.color('marine', ep_site)
pymol.viewing.show('mesh',ct_site)
pymol.viewing.color('tv_blue',ct_site)
pymol.viewing.show('mesh', f267_site)
pymol.viewing.color('violet', f267_site)
pymol.viewing.show('mesh', w336_site)
pymol.viewing.color('teal', w336_site)

# Also show other allosteric sites:
all_allo_sites = ['423', '522', '506', '420']
showing_allo_sites = [i for i in all_allo_sites if (i not in ip_a)]

for sa in showing_allo_sites:
    rn = '/' + pro_name + '///' + sa
    pymol.cmd.select(rn)
    pymol.viewing.label("sele and n. CA", '"%s%s" %(one_letter[resn],resi)')
    center_of_mass.com('sele', object=sa)
    pymol.editing.alter(sa, "vdw=1.4")
    pymol.viewing.show('spheres', sa)
    pymol.viewing.color('salmon', sa)

pymol.cmd.set('label_position','(0,0,10)')
pymol.cmd.set('label_size',16)


# Set view
view = "\
     0.546794116,   -0.278576851,   -0.789562523,\
    -0.586088717,   -0.800804913,   -0.123340204,\
    -0.597923934,    0.530196965,   -0.601146221,\
     0.000000000,    0.000000000, -164.443664551,\
   -59.844493866,  -14.433345795,  -18.295627594,\
   131.734619141,  197.152725220,  -20.000000000 "

pymol.cmd.set_view(view)


# END


# ********** Custimize pathway cylinder ********** #
## For people who want to draw a pathway with fixed cylinder radius,
# or change colors of pathway cylinder. Please use codes provided below:

# # Please enter all residues of a pathway in order.
# pathway = ['420','419','380','381','383','498','499','496']
#
# # Main
# cylinder_pairs = []
# for i in range(0 ,len(using_path) -1):
#     start = using_path[i]
#     end = using_path[i+1]
#     cylinder_pairs.append([start,end])
#
# for i in range(0,len(cylinder_pairs)):
#     start_resi = cylinder_pairs[i][0]
#     end_resi = cylinder_pairs[i][1]
#     name = start_resi + 'to' + end_resi
#     # draw cylinder
#     xyz1 = []
#     pymol.cmd.iterate_state(1, start_resi, 'xyz1.append([x,y,z])', space=locals(), atomic=0)
#     xyz1 = xyz1[0]
#     x1 = xyz1[0]
#     y1 = xyz1[1]
#     z1 = xyz1[2]
#
#     r1,g1,b1 = 0.2,0.6,0.2 # define start color
#     xyz2 = []
#     pymol.cmd.iterate_state(1, end_resi, 'xyz2.append([x,y,z])', space=locals(), atomic=0)
#     xyz2 = xyz2[0]
#     print(xyz2)
#     x2 = xyz2[0]
#     y2 = xyz2[1]
#     z2 = xyz2[2]
#
#     r2,g2,b2 = 0.2,0.6,0.2 # define end color
#     radius = 0.5
#     pymol.cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], name )
#
