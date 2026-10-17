# Core03: two coherent BO electronic channels

This is a minimal extension of the SAME Core02 runtime: generic Hamiltonian
propagation/energy and state-control hooks, state-owned pulse parameters and
work, multicomponent shape, checkpoint fields. Original H2/H3 adapters remain
unchanged. No molecular identity is hardcoded in the runtime.

H2+ inputs must be computed from PySCF1e Coulomb integrals, generalized
Hermitian eigenvectors, AO inversion parity, Lz-action sigma identification
(rejecting pi), actual electronic dipole -r at the nuclear midpoint, and
cross-geometry overlap phase continuation. Nuclear repulsion1/R is included
in both curves. Basis checks are aug-cc-pVTZ vs aug-cc-pVQZ. No SCF/FCI
qualification is inherited from earlier failed electronic solvers.

Nuclear states have shape(2,N) with sum_c integral|chi_c|²dR=1. Local real
Hermitian2x2 exponentials preserve coherence; DST-I gives positive-interval
Dirichlet nuclear boundaries. Midpoint Strang and independently integrated
power <partial H/partial t> are checked by timestep/grid sensitivity. The
prescribed classical finite sin² pulse has an actual physical field amplitude,
intensity, frequency, duration and polarization; comparisons hold frequency
fixed, not retuned for each basis. State-owned controls prevent fork leakage.

The supplied external reference is
https://arxiv.org/abs/1705.03608 . Its table1 energy gaps are Ry, converted to
Eh by division by2. Only attributed numerical facts are included; no full
article republication. PySCF/libcint integral definitions are documented at
https://pyscf.org/_modules/pyscf/gto/moleintor.html . PySCF is an installed
dependency, not redistributed source; original QML MIT notices remain.

This model freezes molecular rotation and COM and includes only two BO
electronic channels. It is NOT intrinsic nonadiabatic dynamics, ionization,
many-electron spin coupling, free3D dynamics or a physical discovery. Finite
basis/incomplete channels/interpolation can limit actual physical accuracy.
A fixed-R analytic Rabi test verifies an exponent, not real dipole accuracy.

`QuantumWorld(system).prepare/apply/evolve/observe/fork/compare/checkpoint/
export_evidence` is the same API. `requests.json` supplies runtime conditions.
The view extends Core02 components and accepts electric-pulse JSON. It is
assembled ONLY after actual electronic/dynamic data exist, never from mock
surfaces or invented dipoles. Node parity is distinct from visual browser QA.

External write/run authority is NOT granted by the scientific design alone.
The exact public file hashes and one-run scope are prepared locally for one
Human approval. No Core03 GitHub push/dispatch or real PySCF calculation has
yet occurred. Mac/Hetzner high-load research NO-GO remains unchanged.
