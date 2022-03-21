""" --- Notes ---

This script is for summarising AutoDock4 and PLIP information of docking poses.

    --Data--
     Materials:
    PLIP results

     Require inputs:
     1. Specify paths at the beginning of script.

     Output:
     1. Summary txt for each ligand+allosteric combination

"""

import os

# Set paths
path_plip = "/outputs/PLIP/" #Please change accordingly
os.chdir(path_plip)
all_files = os.listdir()

all_fd_names = []
all_fd_paths = []
for file in all_files:
    if os.path.isdir(file):
        all_fd_names.append(file)
        all_fd_paths.append(os.getcwd() + '/' + file)
print(all_fd_names)
print(all_fd_paths)

ligand_lists = ['9NO2OA', '10NO2OA', '10NO2LA','PGJ2']
site_lists = ['423', '522', '477',
                      '518', '513', '506', '420']
chain_lists = ['A']
type_lists = ['Hydrophobic Interactions', 'Salt Bridges', 'Hydrogen Bonds']

#TODO 1: for each ligand - for each site -obtain ligand/chain/site/rank/run infos
# for lig, site, chain in zip(ligand_lists, site_lists, chain_lists):
for lig in ligand_lists:
    for site in site_lists:
        for chain in chain_lists:
            #generate general lists for each run & print heading
            pc_fd = []
            pc_ranks = []
            pc_runs = []
            pc_energies = []
            pc_tps = []
            pc_hb = []
            pc_hpho = []
            pc_sb = []
            pc_others = []
            pc_resi_A = []

            pc_n = lig + '_' + chain + '_' + site
            print('Processing: ' + pc_n)
            for fd_n in all_fd_names:
                if fd_n.startswith(pc_n):
                    pc_fd.append(fd_n)

            if bool(pc_fd) is False:
                print('No files available for ' + pc_n)
                continue
            else:
                print('Found ' + str(len(pc_fd)) + ' folders for ' + pc_n)

                for fd in pc_fd:
                    fd_info = fd.replace((pc_n + '_'),'').split('_')
                    #sample fd_info(not split): docking_rank5_run5_energyn3p41_for_interaction
                    rank = fd_info[1].replace('rank','')
                    run = fd_info[2].replace('run','')
                    energy = fd_info[3].replace('m','-').replace('p','.').replace('energy','')
                    rp_path = os.getcwd() + '/' + fd + '/report.txt'
                    hb_n = 0
                    hpho_n = 0
                    sb_n = 0
                    others_n = 0
                    interacting_resi = []
                    tps = []
                    #open report files and generate info
                    with open(rp_path,'r') as f:
                        lines = f.readlines()
                        #check types of interactions
                        for line in lines:
                            if line.startswith('**'):
                                tps.append(line.replace('**','').replace('\n',''))
                            # print(tps)
                        # for each types
                        for tp in tps:
                            recording = False
                            tp_resi = []
                            if ('Hydrophobic Interactions') in str(tp):
                                for line in lines:
                                    # print(recording)
                                    if recording is False:
                                        if str('**' + tp) in str(line):
                                            recording = True
                                    else:
                                        if line.startswith('|') & (('RESNR') not in line):
                                            resi = line.replace(' ', '').replace('/n', '').split(sep='|')[1]
                                            # print(resi)
                                            tp_resi.append(str(resi))
                                        elif not line.strip():
                                            recording = False
                                            break
                                hpho_n = len(tp_resi)
                                interacting_resi.extend(tp_resi)
                                # print(hpho_n)
                            elif ('Salt Bridges') in str(tp):
                                for line in lines:
                                    # print(recording)
                                    if recording is False:
                                        if str('**' + tp) in str(line):
                                            recording = True
                                    else:
                                        if line.startswith('|') & (('RESNR') not in line):
                                            resi = line.replace(' ', '').replace('/n', '').split(sep='|')[1]
                                            # print(resi)
                                            tp_resi.append(str(resi))
                                        elif not line.strip():
                                            recording = False
                                            break
                                sb_n = len(tp_resi)
                                interacting_resi.extend(tp_resi)
                                # print(sb_n)
                            elif ('Hydrogen Bonds') in str(tp):
                                for line in lines:
                                    # print(recording)
                                    if recording is False:
                                        if str('**' + tp) in str(line):
                                            recording = True
                                    else:
                                        if line.startswith('|') & (('RESNR') not in line):
                                            resi = line.replace(' ', '').replace('/n', '').split(sep='|')[1]
                                            # print(resi)
                                            tp_resi.append(str(resi))
                                        elif not line.strip():
                                            recording = False
                                            break
                                hb_n = len(tp_resi)
                                interacting_resi.extend(tp_resi)
                                # print(hb_n)
                            else: #NOTES if there are more than 1 others tp, this section need to be adjusted
                                for line in lines:
                                    # print(recording)
                                    if recording is False:
                                        if str('**' + tp) in str(line):
                                            recording = True
                                    else:
                                        if line.startswith('|') & (('RESNR') not in line):
                                            resi = line.replace(' ', '').replace('/n', '').split(sep='|')[1]
                                            # print(resi)
                                            tp_resi.append(str(resi))
                                        elif not line.strip():
                                            recording = False
                                            break
                                others_n = len(tp_resi)
                                interacting_resi.extend(tp_resi)

                    if site not in interacting_resi:
                        interacting_resi.append(site)

                    inter_A = []
                    for inter in interacting_resi:
                        inter_A.append('A.' + str(inter))
                    input_A = ', '.join(sorted(set(inter_A)))  # set == distinct

                    pc_ranks.append(rank)
                    pc_runs.append(run)
                    pc_energies.append(energy)
                    pc_tps.append(tps)
                    pc_hb.append(hb_n)
                    pc_hpho.append(hpho_n)
                    pc_sb.append(sb_n)
                    pc_others.append(others_n)
                    pc_resi_A.append(input_A)


            #write summary
            with open(os.getcwd() + '/summary_' + pc_n + '.txt', 'w') as wf:
                wf.write('rank; run; energy; types; h_bonds; hydropho_bonds; saltbridges; others; input_A \n')
                for (prk,prn,pce,ptp,phb,phpho,psb,pot,pinpA) in zip(pc_ranks, pc_runs,pc_energies ,pc_tps,
                                                    pc_hb, pc_hpho, pc_sb, pc_others, pc_resi_A):
                    wf.write("{0}; {1}; {2}; {3}; {4}; {5}; {6}; {7}; {8} \n".format(prk,prn,pce,ptp,phb,phpho,psb,pot,pinpA))
                wf.close()

