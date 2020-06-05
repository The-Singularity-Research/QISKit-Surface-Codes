#
# def x_measurement(circuit: qiskit.circuit.QuantumCircuit, qubit: int, cbit: int):
#     """Measure 'qubit' in the X-basis, and store the result in 'cbit'"""
#     circuit.measure = measure  # fix a bug in qiskit.circuit.measure
#     circuit.h(qubit)
#     circuit.measure(qubit, cbit)
#     circuit.h(qubit)
#     return circuit
#
#
# def compute_phi(sigma: Tuple[Tuple[int]], alpha: Tuple[Tuple[int]]) -> List[List[int]]:
#     """compute the list of lists full cyclic form of phi (faces of dessin [sigma, alpha, phi])"""
#     s = Permutation(sigma)
#     a = Permutation(alpha)
#     f = ~(a * s)
#     f = f.full_cyclic_form
#     return f
#

def permlist_to_tuple(perms):
    """
    convert list of lists to tuple of tuples in order to have two level iterables
    that are hashable for the dictionaries used later
    """
    return tuple(tuple(perm) for perm in perms)
#
#
# def surface_code_graph(sigma: Tuple[Tuple[int]], alpha: Tuple[Tuple[int]]) -> \
#         Tuple[Union[nx.MultiGraph, Any], Tuple[Dict[Any, int], Dict[Any, int], Dict[Any, int]]]:
#     """output a graph and a node_info"""
#     f = compute_phi(sigma, alpha)
#     phi = permlist_to_tuple(f)
#
#     surface_graph = nx.MultiGraph()
#
#     node_info = build_node_info(sigma, alpha, phi)
#
#     # Create nodes for each cycle in sigma
#     for cycle in sigma:
#         surface_graph.add_node(cycle, bipartite=1)
#         for node in cycle:
#             surface_graph.add_node(node, bipartite=0)
#             surface_graph.add_edge(cycle, node)
#
#     # Create nodes for each cycle in phi
#     for cycle in phi:
#         surface_graph.add_node(cycle, bipartite=1)
#         for node in cycle:
#             surface_graph.add_edge(cycle, node)
#
#     # Create nodes for each cycle in alpha
#     # then glue the nodes corresponding to a the pairs
#     for pair in alpha:
#         surface_graph.add_node(pair)
#         surface_graph = nx.contracted_nodes(surface_graph, pair[0], pair[1], self_loops=True)
#         surface_graph = nx.contracted_nodes(surface_graph, pair, pair[0], self_loops=True)
#
#     return surface_graph, node_info
#
#
# def build_node_info(sigma, alpha, phi):
#     sigma_dict = dict()
#     count = -1
#     for count, cycle in enumerate(sigma):
#         sigma_dict[cycle] = count
#     phi_dict = dict()
#     for count, cycle in enumerate(phi, start=count + 1):
#         phi_dict[cycle] = count
#     alpha_dict = dict()
#     for count, pair in enumerate(alpha, start=count + 1):
#         alpha_dict[pair] = count
#     return sigma_dict, alpha_dict, phi_dict
#
#
# def surface_code_circuit(sigma: Tuple[Tuple[int]], alpha: Tuple[Tuple[int]]) -> \
#         Tuple[QuantumCircuit, Union[nx.MultiGraph, Any], Tuple[Dict[Any, int], Dict[Any, int], Dict[Any, int]]]:
#     surface_graph, node_info = surface_code_graph(sigma, alpha)
#
#     '''
#     Compute the permutation corresponding to phi and create a
#     'surface code circuit' based on a (multi)graph 'surface_graph'
#     given by sigma, alpha, and phi
#     Create quantum and classical registers based on the number of nodes in G
#     '''
#     f = compute_phi(sigma, alpha)
#     phi = permlist_to_tuple(f)
#
#     qr = QuantumRegister(len(surface_graph.nodes))
#     cr = ClassicalRegister(len(surface_graph.nodes))
#     circ = QuantumCircuit(qr, cr)
#     sigma_dict, phi_dict, alpha_dict = node_info
#     for cycle in sigma:
#         circ.h(sigma_dict[cycle])
#
#     for cycle in phi:
#         circ.h(phi_dict[cycle])
#
#     return circ, surface_graph, node_info
#
#
# def face_syndrome_measure(circuit: qiskit.circuit.QuantumCircuit, surface_graph: nx.MultiGraph,
#                           node_info: Tuple[Dict[Tuple[int], int], Dict[Tuple[int], int], Dict[Tuple[int], int]],
#                           vertex: Tuple[int]):
#     v = vertex
#     sigma_dict, phi_dict, alph_dict = node_info
#     for node in surface_graph.neighbors(v):
#         circuit.cz(phi_dict[v], alph_dict[node])
#
#     circuit.barrier()
#     x_measurement(circuit, phi_dict[v], phi_dict[v])
#     circuit.barrier()
#
#     return circuit, surface_graph, node_info
#
#
# def star_syndrome_measure(circuit: qiskit.circuit.QuantumCircuit, surface_graph: nx.MultiGraph,
#                           node_info: Tuple[dict[Tuple[int]:int]], vertex: Tuple[int]):
#     v = vertex
#     sigma_dict, phi_dict, alph_dict = node_info
#     for node in surface_graph.neighbors(v):
#         circuit.cx(sigma_dict[v], alph_dict[node])
#
#     circuit.barrier()
#     x_measurement(circuit, sigma_dict[v], sigma_dict[v])
#     circuit.barrier()
#
#     return circuit, surface_graph, node_info
