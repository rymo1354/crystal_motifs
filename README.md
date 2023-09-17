# crystal_motifs
- Tools for characterizing and comparing different sub-structure motifs between periodic crystal structures
- Scripts to classify DFT optimized inorganic and organic materials as 3D perovskites according to their octahedral network connectivities
- Handles multiple inorganic A and B-sites, and constructs organic A-sites from periodic structures 

# dependencies
- pymatgen version 2022.0.5 or greater (for pymatgen nearest neighbors finder compatibility)
- scipy (convex hull construction)
- networkx (graph generation)
- pubchempy (querying organic molecules)
- openbabel (smile string representations of organics)

pubchempy and openbabel used to identify molecular perovskite sites, see MAPbI3 test. 

# classify structures as 3D perovksites on HPC resources

Example execution: ```classify_submitter.py -s classify_3D.py -d {DIRECTORY_NAME} -w {FILENAME.json} -o {FILENAME} -n 1 -t 1 -cu```

- ```-s```: classification script to run. ```classify_3D.py``` classifies structures in parallel on the total number of cores available based on the ```-n``` argument.
- ```-d```: directory with structures, typically a directory tree of VASP optimized structures with corresponding VASP runfiles. 
- ```-w```: .json file where the perovskite/non-perovskite classification is written. 
- ```-o```: file to write the standard output of ```classify_3D.py```, 
- ```-n```: number of supercomputing nodes used to classify structures in parallel.
- ```-t```: wall time. For Eagle supercomputing cluster, wall time = 1 submits to the debug queue.  
- ```-cu```: check unconverged, tag that determines whether unconverged VASP runs should be included with the classification (acts as True/False). 

Note: The script ```classify_3D.py``` currently only checks for files named "CONTCAR" to perform the classification on. <br />
Note: The classification dictionary keys are also each structures' reduced formula. If there are multiple structures in the same subdirectory with the same reduced formula, these can be overwritten. 

# example classification file

{"CsGeI3": {"assignment_1": {"is_3d_perovskite": true, "A": ["Cs+"], "B": ["Ge2+"], "X": ["I-"]}, "assignment_2": {"is_3d_perovskite": false, "A": ["Ge2+"], "B": ["Cs+"], "X": ["I-"]}}, <br />
"CsPbBr3": {"assignment_1": {"is_3d_perovskite": true, "A": ["Cs+"], "B": ["Pb2+"], "X": ["Br-"]}, "assignment_2": {"is_3d_perovskite": false, "A": ["Pb2+"], "B": ["Cs+"], "X": ["Br-"]}}, <br />
"CsPbI3": {"assignment_1": {"is_3d_perovskite": true, "A": ["Cs+"], "B": ["Pb2+"], "X": ["I-"]}, "assignment_2": {"is_3d_perovskite": false, "A": ["Pb2+"], "B": ["Cs+"], "X": ["I-"]}}, <br />
"RbPbI3": {"assignment_1": {"is_3d_perovskite": true, "A": ["Rb+"], "B": ["Pb2+"], "X": ["I-"]}, "assignment_2": {"is_3d_perovskite": false, "A": ["Pb2+"], "B": ["Rb+"], "X": ["I-"]}}, <br />
"CsSnBr3": {"assignment_1": {"is_3d_perovskite": true, "A": ["Cs+"], "B": ["Sn2+"], "X": ["Br-"]}, "assignment_2": {"is_3d_perovskite": false, "A": ["Sn2+"], "B": ["Cs+"], "X": ["Br-"]}}, <br />
"CsSnI3": {"assignment_1": {"is_3d_perovskite": true, "A": ["Cs+"], "B": ["Sn2+"], "X": ["I-"]}, "assignment_2": {"is_3d_perovskite": false, "A": ["Sn2+"], "B": ["Cs+"], "X": ["I-"]}}} <br />

# example output file 

Running analysis... <br />
CsPbI3: Not all Cs+ B-site nodes are octahedrons <br />
CsPbBr3: Not all Cs+ B-site nodes are octahedrons <br />
CsSnBr3: Not all Cs+ B-site nodes are octahedrons <br />
RbPbI3: Not all Rb+ B-site nodes are octahedrons <br />
CsSnI3: Not all Cs+ B-site nodes are octahedrons <br />
CsGeI3: Not all Cs+ B-site nodes are octahedrons <br />
{'CsGeI3': {'assignment_1': {'is_3d_perovskite': True, 'A': ['Cs+'], 'B': ['Ge2+'], 'X': ['I-']}, 'assignment_2': {'is_3d_perovskite': False, 'A': ['Ge2+'], 'B': ['Cs+'], 'X': ['I-']}}, 'CsPbBr3': {'assignment_1': {'is_3d_perovskite': True, 'A': ['Cs+'], 'B': ['Pb2+'], 'X': ['Br-']}, 'assignment_2': {'is_3d_perovskite': False, 'A': ['Pb2+'], 'B': ['Cs+'], 'X': ['Br-']}}, 'CsPbI3': {'assignment_1': {'is_3d_perovskite': True, 'A': ['Cs+'], 'B': ['Pb2+'], 'X': ['I-']}, 'assignment_2': {'is_3d_perovskite': False, 'A': ['Pb2+'], 'B': ['Cs+'], 'X': ['I-']}}, 'RbPbI3': {'assignment_1': {'is_3d_perovskite': True, 'A': ['Rb+'], 'B': ['Pb2+'], 'X': ['I-']}, 'assignment_2': {'is_3d_perovskite': False, 'A': ['Pb2+'], 'B': ['Rb+'], 'X': ['I-']}}, 'CsSnBr3': {'assignment_1': {'is_3d_perovskite': True, 'A': ['Cs+'], 'B': ['Sn2+'], 'X': ['Br-']}, 'assignment_2': {'is_3d_perovskite': False, 'A': ['Sn2+'], 'B': ['Cs+'], 'X': ['Br-']}}, 'CsSnI3': {'assignment_1': {'is_3d_perovskite': True, 'A': ['Cs+'], 'B': ['Sn2+'], 'X': ['I-']}, 'assignment_2': {'is_3d_perovskite': False, 'A': ['Sn2+'], 'B': ['Cs+'], 'X': ['I-']}}} <br />
Writing data... <br />
Elapsed time: 19.94055724143982 <br />
Process time: 5.400594236 <br />


