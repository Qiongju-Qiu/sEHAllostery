# sEHAllostery
Allosteric regulation of the soluble epoxide hydrolase by nitro fatty acids by a combined experimental and computational approach

### Introduction 
***
This repos provides raw data and analytic scripts for in silico analysis of corresponding paper.
For more details please see NOTES in each script.

### Details
***
To reproduce in silico analysis, you will need: （brackets indicate locations of corresponding files)

| Studies | provided materials | additional requirements |
| :----------: | :-----------:  | :-----------: |
| Docking | protein & ligands structures (general); <br>Docking_histidine_pairfit.py (python); <br>Docking_pose.py (python) | Python 2 or 3, PyMOL; <br>Scripts and other required software for AutoDock4 Covalent docking (require Python 2): https://autodock.scripps.edu/resources/covalent-docking/|
| Docking results processing | Docking results (outputs/docking/for_PLIP); <br>Docking_record_interaction.py (python)  | Python 3; <br>PLIP (https://github.com/pharmai/plip) |
| Allosteric signalling analysis | Allosteric signalling anaysis raw data (raw_data/AlloSigMA/allosteric_signalling_analysis) | - |
| Allosteric signalling map (ASM) analysis | ASM raw_data (raw_data/AlloSigMA/ASM); <br>ASM_json.py and ASM_signleplot.R (python) | Python 3, R|
| OHM analysis | OHM raw_data (raw_data/OHM); <br>OHM_pathway.py and center_of_mass.py (python) | Python 3, PyMOL |


