# Surface Codes and Error Correction

A good introductory video on closely related work at IBM can be found here: 
- [A graph-based formalism for surface codes and twists](https://www.youtube.com/watch?v=Ca85qdptceQ)

For an introduction on how to use this library see the following notebook:

- [Surface Code Tutorial II](https://github.com/The-Singularity-Research/QISKit-Surface-Codes/blob/master/surface_code_tutorial_2.ipynb)

## Current Capabilities 

This package implements general "*surface codes*" which 
generalize A. Kitaev's toric code for quantum error 
correction. In particular, any code specified by a 
bipartite graph cellularly embedded in a compact Riemann
surface can be constructed using this package. Obviously error
correction is necessary for fault tolerant quantum computing
and surface codes are one of the best studied, most
robust, and most easily implemented types of error 
correction. Furthermore, they only require nearest neighbor
qubit interaction (gates) for implementation. The surface codes
implemented in this package include the toric code as well as
the hyperbolic surface codes described in:
- [Constructions and Noise Threshold of Hyperbolic
Surface Codes](https://arxiv.org/pdf/1506.04029.pdf)
- [Homological Quantum Codes
Beyond the Toric Code](https://arxiv.org/pdf/1802.01520.pdf)

These codes provide a higher error rate threshold than the 
standard toric codes, which have an error threshold of around 
1%. The codes can be transpiled in QISKit and run on any IBM 
hardware backend, or they can be run on specialized hardware 
for hyperbolic surface code such as the hardware constructed 
with around 150 qubits in:
- [Hyperbolic Lattices in Circuit Quantum Electrodynamics](https://arxiv.org/pdf/1802.09549.pdf)

and used in the simulation of hyperbolic spaces in 
- [Quantum Simulation of Hyperbolic Space with Circuit Quantum Electrodynamics:
From Graphs to Geometry](https://arxiv.org/pdf/1910.12318.pdf)

With error rates on current hardware now reaching the 0.1% range
and error thesholds well above 1% using surface codes, 
the implementation of error correction using surface codes is both
timely and necessary for near term applications and will provide
fault tolerance on most current hardware. The primary reference 
for implementing error correction protocols and devloping 
applications to condensed matter physics and spin-glass models is: 
- [Quantum Computation with Topological Codes: from qubit to topological fault-tolerance](https://arxiv.org/pdf/1504.01444.pdf)

## Future Applications 
Future applications will involve: 
- modeling zeros of L-functions given by partition functions 
of lattice (or graph) Ising type models in order to understand 
the behavior of prime numbers and RSA public key cryptography 
and elliptic curve cryptography, 
- understanding the Yang-Mills mass gap problem through 
modeling lattice gauge theories using surface codes, 
- modeling arbitrage and cryptocurrency economies using lattice 
gauge theory, 
- training machine learning models on simulations of these 
models to predict arbitrage vulnerabilities and potential 
cryptocurrency instabilities
- training machine learning models on exotic quantum hardware 
architecture to learn and improve error correction protocols 
and logical gate implementations on non-standard surface codes 
that go beyond simple grid qubit layouts on standard quantum 
hardware, 
- modeling quantum gravity(ies) and the AdS/CFT correspondecne 
of J. Maldecena
- modeling quantum dynamical systems to better understand
quantum chaos, quantum complexity, and phase transitions in
order to better understand molecular dynamics
- using simulations of molecular dynamics to to train 
machine learning models for drug discovery and materials
discovery
- using molecular dynamics to understand complex quantum 
phenomena such as protein folding, molecular biology processes,
genetics for CRISPR, and synthetic biology

## Companion Work
This package will be a part of a larger whole on error correction. 
The error correction software will also implement graph states
and basic applications of graph states such as:

- encoding arbitrary stabilizer codes as graph states
- Measurement Based Quantum Computing (MBQC)
- Bipartite graph states and quantum cryptography
- blind quantum computation
- modeling entanglement entropy of multipartite systems
- computing stabilizer generators of arbitrary stabilizer codes
- weighted graph states and generalized graph states with
arbitrary (variational) controlled-U3 gates as entangling 
gates for quantum machine learning on weighted graphs and 
graph completion problems in network analysis and 
knowledge graphs

Currently some basic constructors of graph states based on 
arbitrary input graphs, as well as bipartite input graphs for
MBQC, quantum cryptography, and blind quantum computation have
been implemented, along with some basic measurement protocols.
This should be ready to use in the very near future as well.

IBM has already implemented simple repetition codes and will 
likely implement some of the standard examples of stabilizer 
codes such as Calderbank-Shor-Steane codes (CSS-codes), so we 
will likely not focus on implementing those in their standard
form. However, developing algorithms which convert these codes
and arbitrary stabilizer codes efficiently into graph states 
or surface codes will likely be useful and is something we are 
likely to focus on in the future.

For supplementary material explaining some of the theoretical 
aspects of surface codes, see the following [notebook](https://github.com/The-Singularity-Research/Surface-Codes)

