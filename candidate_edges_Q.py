from qiskit import Aer, assemble, QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import Statevector
import numpy as np
from qiskit.providers.aer import AerSimulator
from qiskit.quantum_info.operators import Operator
from qiskit.circuit.library.phase_oracle import PhaseOracle
from qiskit.circuit.library import MCMT, UnitaryGate, HGate
from qiskit.result import QuasiDistribution
import IPython
import random
import time

# works only on classical input
# input and expected_output is written starting from the leftmost qubit register of the gate (reverse of the circuit qubit order)
# input is an array of input qubits values (0/1) for the gate (reverse of the circuit qubit order)
# expected_output is a string composed of expected values of qubits listed in expected_output_qubits
# expected_output_qubits is based on the circuit qubit order (see example)
def test_gate(gate, input, expected_output, expected_output_qubits):
    # example: test_gate(gate, [b_1, b_0, a_3, a_2, a_1, a_0], '11', [b_1_qubit_index, b_0_qubit_index])
    # -> test_gate(gate, [0, 0, 1, 1, 0, 1], '11', [4, 5]) # notice the expected_output_qubits indexes
    
    q_r = QuantumRegister(gate.num_qubits, 'q')
    cl_r = ClassicalRegister(len(expected_output_qubits))
    
    circuit = QuantumCircuit(q_r, cl_r)

    for i in range(len(input)):
        circuit.initialize(input[len(input)-1 - i], i)

    circuit.append(gate, q_r[:])
    for i in range(len(expected_output_qubits)):
        circuit.measure(expected_output_qubits[i], i)

    #print(circuit)

    simulator = Aer.get_backend('qasm_simulator')
    t_circ = transpile(circuit, simulator)
    result = simulator.run(t_circ).result()
    counts = result.get_counts()

    # IPython.display.display(plot_histogram(counts))

    # TODO: extend this part of the code to work on non-classical input

    print(counts)
    
    if  (expected_output in counts) and (counts[expected_output] == 1024):
        print("Test successful!")
    else:
        print("Test failed!")


# |d>|c>|b>|a> -> |y>|c>|b>|a>
def get_candidate_edges_gate(point_index_size, query_func, grid_width, direction, number_of_true_points):
    # grid_width -> number of nodes in the first row of the grid
    # direction = 0 -> right-down
    # direction = 1 -> right
    # direction = 2 -> left-down
    
    if len(query_func) != 2 ** point_index_size:
        print("Length of the query function does not match the size of a point index!")
        return
    
    a_size = point_index_size
    b_size = 2

    # U_b
    
    U_b_mat = np.zeros((2 ** (b_size + a_size), 2 ** (b_size + a_size)), dtype=int)

    # |b>|a> -> |b XOR f(a)>|a>
    for a_val in range(2 ** a_size):
        # a_val is a value of a in dec form
        out_val = query_func[a_val]
        for b_val in range(b_size ** 2):
            mat_column = b_val * (2**a_size) + a_val  # b_val shifted left by `a` size and then ORed with `a` gives `b` concatenated (in binary) with `a`
            mat_row = (b_val ^ out_val) * (2**a_size) + a_val
            U_b_mat[mat_row, mat_column] = 1

    U_b_gate = UnitaryGate(U_b_mat, label='U_b')

    # U_c
    
    c_size = b_size

    U_c_mat = np.zeros((2 ** (b_size + a_size), 2 ** (b_size + a_size)), dtype=int)

    # |c>|a> -> |c XOR f(a)>|a>
    row = 0  # index of a row (%2 = 0 -> number of points = width; %2 = 1 -> number of points = width + 1) 
    row_node = 0  # index of a row node
    for a_val in range(2 ** a_size):
        # a_val is a value of `a` in dec form
        
        if direction == 0:
            offset = grid_width + 1
            if ((row % 2) == 1) and (row_node + 1 == grid_width + 1):  # if there's no right-down edge from this node
                offset = 0
        elif direction == 1:
            offset = 1
            if (row_node + 1) == (grid_width + (row % 2)):  # if the node is on the right boundary
                offset = 0
        elif direction == 2:
            offset = grid_width
            if ((row % 2) == 1) and (row_node + 1 == 1):  # if there's no left-down edge from this node
                offset = 0

        if ((a_val + offset) >= 2 ** a_size) or ((a_val + offset) >= number_of_true_points):
            offset = 0
        
        out_val = query_func[a_val + offset]
        
        for c_val in range(c_size ** 2):
            mat_column = c_val * (2**a_size) + a_val  # c_val shifted left by `a` size and then ORed with `a` gives `c` concatenated (in binary) with `a`
            mat_row = (c_val ^ out_val) * (2**a_size) + a_val
            U_c_mat[mat_row, mat_column] = 1

        row_node += 1
        if row % 2 == 0:
            if row_node == grid_width:
                row += 1
                row_node = 0
        else:
            if row_node == (grid_width + 1):
                row += 1
                row_node = 0
        

    U_c_gate = UnitaryGate(U_c_mat, label='U_c')

    d_size = 1
    
    U_d_mat = np.zeros((2 ** (d_size + c_size + b_size), 2 ** (d_size + c_size + b_size)), dtype=int)

    # |d>|c>|b> -> |y>|c>|b>
    for b_val in range(2 ** b_size):
        for c_val in range(2 ** c_size):
            for d_val in range(2 ** d_size):
                mat_column_0 = d_val * (2 ** (b_size + c_size)) + (b_val * (2**c_size) + c_val)
                out_val = (1 if abs(b_val - c_val) == 2 else 0)
                mat_row = (d_val ^ out_val) * (2 ** (b_size + c_size)) + (b_val * (2**c_size) + c_val)
                U_d_mat[mat_row, mat_column_0] = 1

    U_d_gate = UnitaryGate(U_d_mat, label='U_d')

    U_b_inverse_gate = UnitaryGate(U_b_mat.T, label='U_b_I')
    U_c_inverse_gate = UnitaryGate(U_c_mat.T, label='U_c_I')

    a_r = QuantumRegister(a_size, 'a')
    b_r = QuantumRegister(b_size, 'b')
    c_r = QuantumRegister(c_size, 'c')
    d_r = QuantumRegister(d_size, 'd')
    
    circuit = QuantumCircuit(a_r, b_r, c_r, d_r)

    circuit.append(U_b_gate, a_r[:] + b_r[:])
    circuit.append(U_c_gate, a_r[:] + c_r[:])
    circuit.append(U_d_gate, b_r[:] + c_r[:] + d_r[:])
    circuit.append(U_c_inverse_gate, a_r[:] + c_r[:])
    circuit.append(U_b_inverse_gate, a_r[:] + b_r[:])
   
    # print(circuit)

    candidate_edges_gate = circuit.to_gate(label='Z_f')

    return candidate_edges_gate

def base_state(number_of_qubits, state_index):
    column_vector = np.zeros((2**number_of_qubits, 1), dtype=complex)
    column_vector[state_index][0] = 1
    return column_vector

def get_state_index(number_of_qubits, state):
    state_index = 0
    exp = 2**(number_of_qubits - 1)
    for k in state:
        state_index += k * exp
        exp /= 2
    return state_index

def grover_iteration(work_qubits, all_qubits, circuit, candidate_edges_gate):
    circuit.append(candidate_edges_gate, all_qubits)
    
    circuit.h(work_qubits)
    
    # Z_OR operator
    zero_state_ket = np.zeros((2**len(work_qubits), 1))
    zero_state_ket[0, 0] = 1 
    zero_state_bra = zero_state_ket.conjugate().transpose()
    Z_OR_mat = 2 * np.dot(zero_state_ket, zero_state_bra) - np.identity(2**len(work_qubits)) # Z_OR = 2*|00..0><00..0| - 1
    
    Z_OR = UnitaryGate(Z_OR_mat, label='Z_OR')
    
    circuit.append(Z_OR, work_qubits)

    circuit.h(work_qubits)

def find_candidate_edges(quadrants_arr, grid_width):
    # random.seed(1)
    random.seed(time.time())
    
    point_index_size = int(np.ceil(np.log2(len(quadrants_arr))))
    
    #  padding with zeros to create additional data for `a` register and substract 1 to get quadrants from 0 to 3
    quadrants_data = np.pad(quadrants_arr, (0, 2**point_index_size - len(quadrants_arr)), 'constant', constant_values=(0, 0)) - 1
    
    candidate_edges_list = []
    for direction in range(3):
        candidate_edges_gate = get_candidate_edges_gate(point_index_size=point_index_size, query_func=quadrants_data, grid_width=grid_width, direction=direction, number_of_true_points=len(quadrants_arr))
        
        a_r = QuantumRegister(point_index_size, 'a')
        b_r = QuantumRegister(2, 'b')
        c_r = QuantumRegister(2, 'c')
        d_r = QuantumRegister(1, 'd')
        cl_a_r = ClassicalRegister(a_r.size, 'a_out')
        circuit = QuantumCircuit(a_r, b_r, c_r, d_r, cl_a_r)
        
        # for phase kickback in Z_f
        circuit.initialize('-', d_r[0])
        
        # # prepare a superposition of work qubits
        circuit.h(a_r[:])
        
        N = 2 ** a_r.size
        # number of solutions
        # s = 4
        # theta = np.arcsin(np.sqrt(s / N))
        # num_of_iterations = np.floor(np.pi / (4 * theta)).astype(int)
        # p = np.sin((2 * num_of_iterations + 1) * theta) ** 2
        # print("P = %12.10f" % (p))
        
        # for unknown s (number of solutions)
        num_of_iterations = random.randint(1, np.floor(np.pi * np.sqrt(N) / 4).astype(int))
        print(f'Number of iterations: {num_of_iterations}')
        
        for i in range(num_of_iterations):
            grover_iteration(a_r[:], a_r[:] + b_r[:] + c_r[:] + d_r[:], circuit, candidate_edges_gate)
            
        for i in range(a_r.size):
            circuit.measure(a_r[i], i)
        
        # print(circuit)
        
        simulator = Aer.get_backend('qasm_simulator')
        t_circ = transpile(circuit, simulator)
        result = simulator.run(t_circ).result()
        counts = result.get_counts()

        # print(f'Direction: {direction}')
        IPython.display.display(plot_histogram(counts))

        reduced_results = []
        max_result = 0
        for result in counts:
            if counts[result] > max_result:
                max_result = counts[result]
        
        for result in counts:
            if counts[result] >= max_result/2:
                reduced_results.append(result)
        #reduced_recults.sort()

        # if over alpha percent of all possible solutions matched, then there's high probability that there was no solution for this direction
        alpha = 0.33

        if len(reduced_results) / 2**a_r.size >= alpha:
            print(f'No solutions for this direction')
        else:
            for node_in_str in reduced_results:
                nodeA = int(node_in_str, 2)
                if direction == 0:
                    nodeB = nodeA + grid_width + 1
                elif direction == 1:
                    nodeB = nodeA + 1
                elif direction == 2:
                    nodeB = nodeA + grid_width
                candidate_edges_list.append([nodeA, nodeB])

    candidate_edges = np.array(candidate_edges_list)
    filter_outside_values = np.logical_and(candidate_edges[:, 0] < len(quadrants_arr), candidate_edges[:, 1] < len(quadrants_arr))

    return candidate_edges[filter_outside_values]














