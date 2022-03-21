""" --- Notes ---

This script is for AutoDock4 covalent docking related to histidine.
It is aim at generate the structures for covalent docking where ligands are attached to corresponding histidine
in a 3D mode.
It was developed using Python 2,in line with scripts provided by AutoDock4.
********* Please run this script in PyMOL. ********

    --Cautions--
    1. Ligands all contain previously-made histidine unit (not necessary).
    2. Before running, change dictionary, ligand indices and residues indices.
    (For definition of ligand indices and residues indices, please refers to readme file of AutoDock4)

    --Data--
     Materials:
     Modified ligand structures

     Require inputs:
     1. path to ligand structures
     2. ligands indices
     3. residues indices

     Output:
     ligand structures with modified coordinates (pdb)

"""

import os
import pymol

#### --- Inputs --- ####
# Input: locate to folder with ligand and protein files
os.chdir("../general")

# define ligand, protein and residues
ligands_list = ['9NO2OA', '10NO2OA', '10NO2LA']

for ligand in ligands_list:
    ligand_path = os.path.join(os.getcwd(),ligand+'_his.mol2')
    if ('9NO2OA') in str(ligand):
        resi_list = ['518', '513', '506', '420', '524']
        carbon = 'C9'
    elif ('10NO2OA') in str(ligand):
        resi_list = ['518', '513', '506', '420', '524']
        carbon = 'C10'
    elif ('10NO2LA_2nd') in str(ligand):
        resi_list = ['513', '506', '524']
        carbon = 'C7'
    elif ('10NO2LA_3rd') in str(ligand):
        resi_list = ['513', '506', '524']
        carbon = 'C6'
    else:
        resi_list = ['513', '506', '524']
        carbon = 'C10'

    for resi in resi_list:
            protein_path = os.path.join(os.getcwd(), 'apoCTDsEH_chain.pdb')

            tg_resi = 'target///HIS`' + resi + '/'
            lig = 'ligand///0`1/'

            lig_rm = ['O32', 'C31', 'C30', 'C29', 'N33', 'C26', 'C25', 'C28', 'N27']


            # Main
            pymol.finish_launching()
            pymol.cmd.load(protein_path, 'target')
            pymol.cmd.load(ligand_path, 'ligand')
            for rm_atom in lig_rm:
                pymol.cmd.edit(lig + rm_atom)
                pymol.cmd.remove_picked()
            pymol.creating.extract('target_resi',tg_resi)
            pymol.cmd.edit(lig + 'N24')
            pymol.cmd.remove_picked()
            pymol.editing.fuse(lig + carbon, 'target_resi///HIS`' + resi + '/NE2')
            pymol.editing.h_add('target_resi')

            at_chain = "'chain=" + '"A"' + "'"
            at_resiv = "'resv=" +  resi  + "'"

            pymol.cmd.alter('target_resi', 'resn="HIS"')
            pymol.cmd.alter('target_resi', eval(at_chain))
            pymol.cmd.alter('target_resi', eval(at_resiv))
            pymol.cmd.alter('target_resi', "type='HETATM'")
            pymol.cmd.alter('target_resi', 'q="1"') #occupancy
            pymol.cmd.alter('target_resi', 'b="0"') #b-factor

            # pymol.fitting.pair_fit(sele)
            pymol.cmd.save(os.getcwd() + '/' + ligand + '_A_his' + resi + '_align.pdb', 'target_resi')

            #CLOSE pymol
            pymol.cmd.reinitialize()

            print('Finished: ' + ligand + ', resi: ' + resi)

