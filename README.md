# GRPF-quantum
A repository for code related to GRPF (Global Root and Pole Finding) algorithm

## How to use
1. **[grpf.py](grpf.py) - starts the program**
2. [analysis_parameters.py](/analysis_parameters.py) - contains all parameters of the analysis, e.g.:
    * the domain shape and size (two domain shapes are available: rectangle and circle, as described in the examples) 
    * the initial step
    * accuracy (Tol)
    * mesh visualization options
3. [fun.py](fun.py) - definition of the function for which roots and poles will be calculated
4. **to run examples**: copy (and overwrite) [analysis_parameters.py](analysis_parameters.py) and [fun.py](fun.py) files from the folder containing the example to the main folder and start GRPF program
 
## Short description of the functions
- [GRPF.py](GRPF.py) - main body of the algorithm  
- [analysis_parameters.py](analysis_parameters.m) - analysis parameters definition
- [fun.py](fun.py) - function definition
- [candidate_edges_Q.py](candidate_edges_Q.py) - generates and simulates quantum circuits that find candidate edges in the first iteration of GRPF
- [rect_dom.py](rect_dom.py) - initial mesh generator for rectangular domain
- [disk_dom.py](disk_dom.py) - initial mesh generator for circular domain
- [vinq.py](vinq.py) - converts the function value into proper quadrant
- [find_next_node.py](find_next_node.py) - finds the next node in the candidate boundary creation process
- [vis.py](vis.py) - mesh visualization