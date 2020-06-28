from typing import Tuple

from qiskit.circuit import QuantumCircuit, QuantumRegister, ClassicalRegister

from SurfaceCodes.surface_code_graph import SurfaceCodeGraph
from SurfaceCodes.utilites import permlist_to_tuple


class SurfaceCodeCircuit(QuantumCircuit):

    def __init__(self, sigma: Tuple[Tuple[int]], alpha: Tuple[Tuple[int]]):
        super().__init__()
        self.sigma = sigma
        self.alpha = alpha

        self.scgraph = SurfaceCodeGraph(self.sigma, self.alpha)

        '''
        Compute the permutation corresponding to phi and create a 
        'surface code circuit' based on a (multi)graph 'surface_code_graph'
        given by sigma, alpha, and phi
        Create quantum and classical registers based on the number of nodes in G
        '''
        # f = self.scgraph.compute_phi()
        self.phi = self.scgraph.phi

        self.qr = QuantumRegister(len(self.scgraph.code_graph.nodes))
        self.cr = ClassicalRegister(len(self.scgraph.code_graph.nodes))
        self.circ = QuantumCircuit(self.qr, self.cr)

        self.node_info = self.scgraph.node_dict 
        self.sigma_dict, self.alpha_dict, self.phi_dict = self.node_info

        for cycle in self.sigma:
            self.circ.h(self.sigma_dict[cycle])

        for cycle in self.phi:
            self.circ.h(self.phi_dict[cycle])

    def x_measurement(self, qubit: int, cbit: int):
        """Measure 'qubit' in the X-basis, and store the result in 'cbit'
        :param qubit, cbit:
        :return None
        """
        # circuit.measure = measure  # fix a bug in qiskit.circuit.measure
        self.circ.h(qubit)
        self.circ.measure(qubit, cbit)
        self.circ.h(qubit)

    def star_syndrome_measure(self, vertex: Tuple[int]):
        """
        Applies CX gates to surrounding qubits of a star then measures star qubit in X-basis
        :param vertex:
        :return:  self.circ, self.scgraph, self.node_info
        """

        for node in self.scgraph.code_graph.neighbors(vertex):
            self.circ.cx(self.sigma_dict[vertex], self.alpha_dict[node])

        self.circ.barrier()
        self.x_measurement(self.sigma_dict[vertex], self.sigma_dict[vertex])
        self.circ.barrier()

        return self.circ, self.scgraph, self.node_info

    def face_syndrome_measure(self, vertex: Tuple[int]):
        """
        Applies CZ gates to surrounding qubits on the boundary of a face then measures face qubit in X-basis
        :param vertex:
        :return:
        """

        for node in self.scgraph.code_graph.neighbors(vertex):
            self.circ.cz(self.phi_dict[vertex], self.alpha_dict[node])

        self.circ.barrier()
        self.x_measurement(self.phi_dict[vertex], self.phi_dict[vertex])
        self.circ.barrier()

        return self.circ, self.scgraph, self.node_info

    def product_Z(self, faces):
        """
        Pauli product Z operator for arbitrary 2-chain boundary
        """

        boundary_nodes = self.scgraph.del_2(faces)
        for node in boundary_nodes:
            self.circ.z(self.alpha_dict[node])

    def product_X(self, stars):
        """
        Pauli product X operator for arbitrary 0-cochain coboundary
        """
        coboundary_nodes = self.scgraph.delta_1(stars)
        for node in coboundary_nodes:
            self.circ.x(self.alpha_dict[node])


    def X_1_chain(self, edges):
        """
        Pauli product X operator for arbitrary 1-cochain given by
        a list of edges
        """
        for edge in edges:
            self.circ.x(self.alpha_dict[edge])


    def Z_1_chain(self, edges):
        """
        Pauli product Z operator for arbitrary 1-chain given by
        a list of edges
        """
        for edge in edges:
            self.circ.z(self.alpha_dict[edge])


    def draw_graph(self, node_type='', layout = ''):
        if layout == 'spring':
            pos=nx.spring_layout(self.scgraph.code_graph)
        if layout == 'spectral':
            pos=nx.spectral_layout(self.scgraph.code_graph)
        if layout == 'planar':
            pos=nx.planar_layout(self.scgraph.code_graph)
        if layout == 'shell':
            pos=nx.shell_layout(self.scgraph.code_graph)
        if layout == 'circular':
            pos=nx.circular_layout(self.scgraph.code_graph)
        if layout == 'spiral':
            pos=nx.spiral_layout(self.scgraph.code_graph)
        if layout == 'random':
            pos=nx.random_layout(self.scgraph.code_graph)
        if node_type == 'cycles':
            self.scgraph.draw('cycles', layout)
        if node_type == 'dict':
            self.scgraph.draw('dict', layout)