""" --- Notes ---

This script is for analysis of autodock4 covalent docking results.
It will extract all top10 docking poses, and generate protein-ligand complex using PyMOL.
The protein-ligand complex will be used as input for interacting residues identifications (with PLIP)
It was developed using Python 3.9.
********* Please run this script in PyMOL. ********

    --Cautions--
    In the output structure of protein-ligand complex,
    some of them shows unexpected C-C bonds between ligand and complex,
    this is caused by setting of PyMOL display and the complex structure is fine.

    Codes for performing PLIP analysis is provided at the end of the script.

    --Data--
     Materials:
    docking results (.dlg)
    protein structure

     Require inputs:
     1. Specify paths at the beginning of script.

     Output:
     1. Ligand structures in folder 'poses'
     2. Protein-ligand structures in folders 'for_PLIP'
     3. PLIP outputs (optional)

"""



import os
import numpy
import pandas as pd
import pymol


# Set paths
# Set path to docking results
path_docking_results = "/outputs/docking/docking_results/"
path_save_ligand_poses = "/outputs/docking/poses/"
path_save_complex = "/outputs/docking/for_PLIP/"
path_protein = "/general/apoCTDsEH.pdb"
os.chdir(path_docking_results)

docking_results_list = []
for r,d,f in os.walk('.'):
    for fname in f:
        if str(fname).endswith('docking.dlg'):
            path = os.path.join(r, fname)
            docking_results_list.append(path)
        else:
            continue


for docking_result in docking_results_list:
    dlg_name = os.path.basename(docking_result).replace('.dlg', '')
    print('Processing: ' + dlg_name)
    dlg_path = docking_result

    if any(resi in str(dlg_name) for resi in ['423', '522', '477']):
        resi_atoms = ['SG', 'CA', 'CB']
    else:
        resi_atoms = ['CA', 'CB', 'CG', 'CD2', 'NE2', 'CE1', 'ND1']

    all_runs = [1,2,3,4,5,6,7,8,9,10]
    energies = []
    all_poses = []
    with open(dlg_path, 'rt') as dlg_f:
        for run in all_runs:
            # print run
            start_line = 'DOCKED: USER    Run = ' + str(run) + '\n'
            recording = False
            pose_info = []
            # print 'start_line: ' + start_line
            with open(dlg_path, 'rt') as dlg_f:
                for line in dlg_f:
                    if recording is False:
                        if line == start_line:
                            pose_info.append(line.strip())
                            recording = True
                    else:
                        if ("ENDMDL") in str(line):
                            pose_info.append(line.strip())
                            recording = False
                        else:
                            # break
                            pose_info.append(line.strip())
            # print pose_info
            best_pose = '\n'.join(pose_info)
            best_pose = best_pose.replace('DOCKED: ', '')
            all_poses.append(best_pose)
            # print best_pose
                    # Record energy and run
            for bp_line in pose_info:
                # print bp_lineEstimated Free Energy of Binding
                if ('Estimated Free Energy of Binding') in bp_line:
                    energies.append(float(bp_line.replace(' kcal/mol  [=(1)+(2)+(3)-(4)]', '').split(' ')[-1].replace('+','')))
        # print energies
        # rank_energies = numpy.array((energies))
        # rank = numpy.argsort(rank_energies)+1
        rank = pd.Series(energies).rank(method='first')
        # print(rank)
        # print len(rank) + len(all_runs) + len(all_poses)
        # print all_poses[1]

        #Outputs: ligand structures
        for p,r,k,e in zip(all_poses, all_runs, rank, energies):
            # #write pdb
            posefin = open(path_save_ligand_poses + dlg_name + '_rank' + str(k) +
                           '_run' + str(r) + '_energy' + str(e) + '.pdb', 'wt')
            posefin.write(p)
            posefin.close()


            # ----pymol main---- #
            pymol.finish_launching()
            pymol.cmd.load(path_save_ligand_poses + dlg_name + '_rank' + str(k) +
                           '_run' + str(r) + '_energy' + str(e) + '.pdb', 'ligand')


            pymol.cmd.alter('ligand', 'segi=""')
            pymol.cmd.alter('ligand', 'resn="lig"')
            pymol.cmd.alter('ligand', 'chain="C"')
            pymol.cmd.alter('ligand', 'resv=0')

            pymol.cmd.remove('hydrogens')
            print('Removed hydrogens for ' + dlg_name)

            for resi_atom in resi_atoms:
                pymol.cmd.edit('ligand//C/lig`0/' + resi_atom)
                pymol.cmd.remove_picked()

            # open corresponding protein file
            pymol.cmd.load(path_protein, 'protein')
            pymol.cmd.select('all')
            pymol.creating.extract('bound_protein', 'sele')

            # save
            pymol.cmd.save(path_save_complex + dlg_name + '_rank' + str(k) +
            '_run' + str(r) + '_energy' + str(e) + '_for_interaction.pdb')
            pymol.cmd.reinitialize()
            #remove atoms
            print('Finished: ' + dlg_name + '_rank' + str(k) + '_run' + str(r))


# Rename complex name so as to perform PLIP in bash shell
os.chdir(path_save_complex)
all_files = os.listdir()
all_pdbs = []
for file in all_files:
    if file.endswith('.pdb'):
        #rename to remove '-' in filename
        os.rename(os.path.join(os.getcwd(),file),os.path.join(os.getcwd(),file.replace('-','m').replace('.0_','_').replace('.pdb','').replace('.','p',1) + '.pdb'))
        all_pdbs.append(file.replace('-','m').replace('.0_','_').replace('.pdb','').replace('.','p',1) + '.pdb')

# # ----perform PLIP---- #
# print(all_pdbs)
# print(len(all_pdbs))
# half = len(all_pdbs)/2
# pdbs_1 = all_pdbs[0:int(half)]
# print(len(pdbs_1))
# pdbs_2 = all_pdbs[int(half):int(len(all_pdbs))]
# print(len(pdbs_2))
#
# # PLIP can perfrom analysis of two pdb file at a time
# for p1,p2 in zip(pdbs_1,pdbs_2):
#     os.system("docker run --rm " + "-v ${PWD}:/results " + "-w /results " + "-u $(id -u ${USER}):$(id -g ${USER}) " +
#               'pharmai/plip:latest -f "' + p1 + '" "' + p2 + '" -yt')
#
