## Docking
### Docking_histidine_pairfit.py 

This script is for AutoDock4 covalent docking related to histidine. <br>
It is aim at generate the structures for covalent docking where ligands are attached to corresponding histidine
in a 3D mode. <br>
It was developed using Python 2,in line with scripts provided by AutoDock4.<br>
<br>
********* Please run this script in PyMOL. ******** <br>
<br>
    **--Cautions--** <br>
    1. Ligands all contain previously-made histidine unit (not necessary).<br>
    2. Before running, change dictionary, ligand indices and residues indices.<br>
    (For definition of ligand indices and residues indices, please refers to readme file of AutoDock4)<br>
    **--Data--**<br>
     **Materials:**<br>
     Modified ligand structures<br>
     <br>
     **Require inputs:**<br>
     1. path to ligand structures<br>
     2. ligands indices<br>
     3. residues indices<br>
     <br>
     **Output:**<br>
     ligand structures with modified coordinates (pdb)<br>

****

 
### Docking_pose.py
This script is for analysis of autodock4 covalent docking results. <br>
It will extract all top10 docking poses, and generate protein-ligand complex using PyMOL. <br>
The protein-ligand complex will be used as input for interacting residues identifications (with PLIP) <br>
It was developed using Python 3.9. <br>
<br>
********* Please run this script in PyMOL. ******** <br>
<br>
    **--Cautions--**<br>
    In the output structure of protein-ligand complex,<br>
    some of them shows unexpected C-C bonds between ligand and complex,<br>
    this is caused by setting of PyMOL display and the complex structure is fine.<br>
    Codes for performing PLIP analysis is provided at the end of the script.<br>
    **--Data--**<br>
    **Materials:**<br>
    docking results (.dlg)<br>
    protein structure<br>
    <br>
     **Require inputs:**<br>
     1. Specify paths at the beginning of script.<br>
     <br>
     **Output:**<br>
     1. Ligand structures in folder 'poses'<br>
     2. Protein-ligand structures in folders 'for_PLIP'<br>
     3. PLIP outputs (optional)<br>
     
   ****  
     
  ### Docking_record_interaction.py
   This script is for summarising AutoDock4 and PLIP information of docking poses. <br>
   <br>
    **--Data--**<br>
     **Materials:**<br>
    PLIP results<br>
    <br>
     **Require inputs:**<br>
     1. Specify paths at the beginning of script. <br>
     <br>
     **Output:**<br>
     1. Summary txt for each ligand+allosteric combination<br>


## Allosteric signalling map (ASM) analysis 
### ASM_json.py
This script is for analysis of Allosteric signalling map (ASM) raw data obtained from AlloSigMA.<br>
It was developed using Python 3.9.<br>
<br>
    **--Cautions--**<br>
     This default setting of this script is to perform analysis at per-residue resolution (one residue at a time)<br>
     For anyone who wants to analyse effects of perturbating two or more residues at the same time,<br>
     please define residues in mt_r_list manually:<br>
     for example, to analyse effects of perturbating 230 and 231 (experiment 1), and 250 and 251 (experiment 2):<br>
     mt_r_list = [['230', '231'], ['250', '251']]<br>
    **--Data--**<br>
     **Materials:**<br>
     ASM raw data<br>
     <br>
     **Require inputs:**<br>
     1. paths for data/outputs;<br>
     2. Specific residues to look at (e.g. residue 383, residue 466);<br>
     3. Specific regions to look at (e.g. Epoxide positioners (383 and 466))<br>
     <br>
     **Output:**<br>
     1.json_all_singular.txt (Effects (delta H) of perturbating one single residue on specific residues;<br>
     2.json_all.txt (Avaerage effects (delta H) of perturbating one signle residue on specific regions.<br>
****

### ASM_singleplot.R
This script is for visulisation of Allosteric Signalling Map (ASM) analysis.<br>
 This script is ONLY suitable for SINGLE residue perturbation, and<br>
to see effects on each functional residue of soluble epoxide hydrolase (sEH).<br>
<br>
For any other uses, please modify accordingly.<br>
 
 **Required inputs:**<br>
  json_all_singular.txt file generated from Python script asm_json.py<br>
   <br>
 **Output:**<br>
   ASM_plot.pdf<br>
<br>

## OHM analysis
### OHM_pathway.py

This script is for analysis and visualisation of OHM pathway data.<br>
It was developed using Python 3.9.<br>
<br>
********* Please run this script in PyMOL. ******** <br>
<br>
    **--Cautions--**<br>
    1. This script provides two types of analysis/visualisation methods:<br>
        # state = 0: display the longest path of top5 pathways;<br>
        # state = 1: display the top1 pathway.<br>
    2. Codes of customize pathway cylinder displays are available at the end of this script.<br>
    <br>
    **--Data--**<br>
     **Materials:**<br>
     OHM data<br>
     center_of_mass.py<br>
     <br>
     **Require inputs:**<br>
     1. Specify path to OHM data.<br>
     2. Specify state of analysis/visualisation.<br>
     3. If customizing pathway cylinder displays, specify all pathway residues.<br>
     <br>
     **Output:**<br>
     PyMOL windows with selected display.<br>
     
