from math import ceil

import matplotlib.pyplot as plt
from qiskit import Aer, assemble, QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import Statevector
import numpy as np
from qiskit.providers.aer import AerSimulator
from qiskit.quantum_info.operators import Operator
from qiskit.circuit.library.phase_oracle import PhaseOracle
from qiskit.circuit.library import MCMT, UnitaryGate, HGate, ZGate, XGate
from qiskit.result import QuasiDistribution
import IPython
import random
import time

from scipy.ndimage import shift

def gate_D(num_b, num_c):
    U_d_mat = np.zeros((2 ** (1 + num_c + num_b), 2 ** (1 + num_c + num_b)), dtype=int)

    b_size = num_b
    c_size = num_c
    d_size = 1
    
    # |d>|c>|b> -> |y>|c>|b>
    for b_val in range(2 ** b_size):
        for c_val in range(2 ** c_size):
            for d_val in range(2 ** d_size):
                mat_column_0 = d_val * (2 ** (b_size + c_size)) + (b_val * (2**c_size) + c_val)
                out_val = (1 if abs(b_val - c_val) == 2 else 0)
                mat_row = (d_val ^ out_val) * (2 ** (b_size + c_size)) + (b_val * (2**c_size) + c_val)
                U_d_mat[mat_row, mat_column_0] = 1
    U_d_gate = UnitaryGate(U_d_mat, label="  $U_d$  ")
    return U_d_gate

def gate_B(num_a, num_b, query_func):
    qc = QuantumCircuit(num_a + num_b)
    x = np.array(list(enumerate(query_func)))
    for n in range(num_b):
        for target in x[:, 0][x[:, 1] & (1 << n) != 0]:
            rev_target = bin(target)[2:].zfill(num_a)[::-1]
            zero_inds = [ind for ind in range(num_a) if rev_target.startswith("0", ind)]

            if zero_inds:
                qc.x(zero_inds)
            qc.compose(MCMT(XGate(), num_a, 1), list(range(num_a)) + [num_a + n], inplace=True)
            if zero_inds:
                qc.x(zero_inds)
    qc.name = "  $U_b$   "
    return qc.to_gate()


def gate_C(num_a, num_c, query_func, direction, grid_width):
    query_func_temp = np.empty_like(query_func)
    index_last_column = len(query_func) - grid_width - (grid_width + 1) % 2
    if direction == 0:
        query_func_temp[:index_last_column - (grid_width % 2)] = query_func[grid_width + 1:]
        query_func_temp[(2 * grid_width)::(2 * grid_width) + 1] = query_func[(2 * grid_width)::(2 * grid_width) + 1]
        query_func_temp[index_last_column:] = query_func[index_last_column:]
    elif direction == 1:
        query_func_temp[:-1] = query_func[1:]
        query_func_temp[grid_width - 1::(2 * grid_width) + 1] = query_func[grid_width - 1::(2 * grid_width) + 1]
        query_func_temp[2 * grid_width::(2 * grid_width) + 1] = query_func[2 * grid_width::(2 * grid_width) + 1]
    elif direction == 2:
        query_func_temp[:index_last_column + (grid_width + 1) % 2] = query_func[grid_width:]
        query_func_temp[grid_width::(2 * grid_width) + 1] = query_func[grid_width::(2 * grid_width) + 1]
        query_func_temp[index_last_column:] = query_func[index_last_column:]

    qc = QuantumCircuit(num_a + num_c)
    qc.append(gate_B(num_a, num_c, query_func_temp), list(range(num_a + num_c)))
    qc.name = "  $U_c$   "
    return qc.to_gate()


# |d>|c>|b>|a> -> |y>|c>|b>|a>
def get_candidate_edges_gate(point_index_size, query_func, grid_width, direction, number_of_true_points):
    # grid_width -> number of nodes in the first row of the grid
    # direction = 0 -> right-down
    # direction = 1 -> right
    # direction = 2 -> left-down

    a_size = point_index_size
    b_size = 2
    c_size = b_size

    gateb = gate_B(a_size, b_size, query_func)
    gatec = gate_C(a_size, c_size, query_func, direction, grid_width)
    gated = gate_D(b_size, c_size)

    a_r = QuantumRegister(a_size, 'a')
    b_r = QuantumRegister(b_size, 'b')
    c_r = QuantumRegister(c_size, 'c')
    d_r = QuantumRegister(1, 'd')

    circuit = QuantumCircuit(a_r, b_r, c_r, d_r)

    circuit.append(gateb, a_r[:] + b_r[:])
    circuit.append(gatec, a_r[:] + c_r[:])
    circuit.append(gated, b_r[:] + c_r[:] + d_r[:])

    gatec_inverse = gatec.inverse()
    gatec_inverse.name = "  $U_c^{-1}$ "

    gateb_inverse = gateb.inverse()
    gateb_inverse.name = "  $U_b^{-1}$ "
    
    circuit.append(gatec_inverse, a_r[:] + c_r[:])
    circuit.append(gateb_inverse, a_r[:] + b_r[:])

    # IPython.display.display(circuit.draw(output='mpl',  style='bw'))

    candidate_edges_gate = circuit.to_gate(label='Z_f')

    return candidate_edges_gate


def base_state(number_of_qubits, state_index):
    column_vector = np.zeros((2 ** number_of_qubits, 1), dtype=complex)
    column_vector[state_index][0] = 1
    return column_vector


def get_state_index(number_of_qubits, state):
    state_index = 0
    exp = 2 ** (number_of_qubits - 1)
    for k in state:
        state_index += k * exp
        exp /= 2
    return state_index


def grover_iteration(work_qubits, all_qubits, circuit, candidate_edges_gate):
    circuit.append(candidate_edges_gate, all_qubits)

    circuit.h(work_qubits)

    # Z_OR operator
    zero_state_ket = np.zeros((2 ** len(work_qubits), 1))
    zero_state_ket[0, 0] = 1
    zero_state_bra = zero_state_ket.conjugate().transpose()
    Z_OR_mat = 2 * np.dot(zero_state_ket, zero_state_bra) - np.identity(
        2 ** len(work_qubits))  # Z_OR = 2*|00..0><00..0| - 1

    Z_OR = UnitaryGate(Z_OR_mat, label='Z_OR')

    circuit.append(Z_OR, work_qubits)

    circuit.h(work_qubits)


def find_candidate_edges(quadrants_arr, grid_width):
    random.seed(time.time())

    point_index_size = int(np.ceil(np.log2(len(quadrants_arr))))  # number of bits to fit a point with largest index
    
    np.set_printoptions(linewidth=(2**point_index_size)*3)

    #  padding with zeros to create additional data for `a` register and substract 1 to get quadrants from 0 to 3
    quadrants_data = quadrants_arr - 1

    candidate_edges_list = []
    for direction in range(3):
        print(f'Direction: {direction}')
        candidate_edges_gate = get_candidate_edges_gate(point_index_size=point_index_size, query_func=quadrants_data,
                                                        grid_width=grid_width, direction=direction,
                                                        number_of_true_points=len(quadrants_arr))

        a_r = QuantumRegister(point_index_size, 'a')
        b_r = QuantumRegister(2, 'b')
        c_r = QuantumRegister(2, 'c')
        d_r = QuantumRegister(1, 'd')
        cl_a_r = ClassicalRegister(a_r.size, 'a_out')
        circuit = QuantumCircuit(a_r, b_r, c_r, d_r, cl_a_r)

        # prepare a superposition of work qubits
        circuit.h(a_r[:])

        N = 2 ** a_r.size
        
        # s = 4  # number of solutions for a given direction (if known)
        # theta = np.arcsin(np.sqrt(s / N))
        # num_of_iterations = np.floor(np.pi / (4 * theta)).astype(int)
        # p = np.sin((2 * num_of_iterations + 1) * theta) ** 2  # probability of finding the solution in one iteration
        # print("P = %12.10f" % (p))

        # for unknown s (number of solutions)
        num_of_iterations = random.randint(1, np.floor(np.pi * np.sqrt(N) / 4).astype(int))

        print(f'Number of iterations: {num_of_iterations}')

        # creating the iterations of grover algorithm
        for i in range(num_of_iterations):
            grover_iteration(a_r[:], a_r[:] + b_r[:] + c_r[:] + d_r[:], circuit, candidate_edges_gate)

        for i in range(a_r.size):
            circuit.measure(a_r[i], i)

        # IPython.display.display(circuit.draw(output='mpl',  style='bw'))
        
        # simulating the circuit 1024 times
        simulator = Aer.get_backend('qasm_simulator')
        t_circ = transpile(circuit, simulator)
        result = simulator.run(t_circ).result()
        counts = result.get_counts()

        if direction == 0:
            text = 'α'
        if direction == 1:
            text = 'β'
        if direction == 2:
            text = 'γ'

        fig = plot_histogram(counts)
        plt.title(f'Direction {text}')
        plt.tick_params(axis='x', which='major', labelsize=3)
        # plt.savefig(f'histogerm_direction_{direction}.eps', format='eps', dpi=1200)
        plt.show()

        reduced_results = []
        max_result = 0
        for result in counts:
            if counts[result] > max_result:
                max_result = counts[result]

        for result in counts:
            if counts[result] >= max_result / 2:
                reduced_results.append(result)

        # if over alpha of all possible solutions matched, then there's high probability that there was no solution for this direction
        alpha = 0.33

        if len(reduced_results) / 2 ** a_r.size >= alpha:
            print(f'No solutions for this direction')
        else:
            # selecting the second point based on the direction of searching
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
    filter_outside_values = np.logical_and(candidate_edges[:, 0] < len(quadrants_arr),
                                           candidate_edges[:, 1] < len(quadrants_arr))

    return candidate_edges[filter_outside_values]
    
