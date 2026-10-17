# Core02: cross-system grid-wavefunction runtime

`runtime.py` has no H2/H3 identity branch. `PhysicalSystem` adapters supply
particles, masses, coordinates, Hamiltonian, boundary, controls, initial states
and model-specific observables. The same symmetric split-operator algorithm
uses either a2D FFT (H3, preserving the mixed kinetic term) or a1D orthonormal
DST-I (H2, reduced radial u(R), Dirichlet). No arbitrary density-matrix or
electronic real-time support is claimed.

Example: `world=QuantumWorld(H2System()); world.allowed();
world.prepare(input_json); world.apply(control_json); world.step(.25);
branch=world.fork(); world.compare(branch); world.export_evidence(path)`.
H3System can replace H2System without changing that call flow. Cross-model
subtraction is rejected. `evolve`, checkpoint and hash-checked restore work on
both. Preparation and controls are model-defined; ideal phase imprints record
work and do not claim a particular laboratory implementation.

H2 is independent of H3/LSTH. Its82 published numerical BO facts come from
Pachucki PRA82,032509(2010), Table4:
https://arxiv.org/html/1007.0322 . Full article text is NOT redistributed.
Cubic interpolation is a model approximation; PCHIP sensitivity, radial-grid
and box changes are reported separately. Published numerical precision does
NOT transfer to interpolated spectra. A separate reference (Dickenson etal,
PRL110,193601(2013), Table2 https://arxiv.org/html/1304.6390) gives the BO
v01 gap4163.40350 cm^-1 and measured gap4161.16632 cm^-1. Missing finite-mass,
relativistic/QED corrections and interpolation error must not be conflated.

The original MIT-licensed QML/LSTH source and original Core01 files remain
unchanged. Only an archived3-step state is used for regression; no original
three full48.38fs experiments are repeated. This branch uses one manual
standard free Ubuntu job; science tree is limited1CPU/1thread/1GiB/no swap/
no network. Container build is preparation on the standard runner, not a
separately hard-limited1GiB host. No artifacts/cache/Packages or paid resources.

The switchable standalone view reuses Core01 FFT/visual styles. Its model
adapters choose preparation, controls, observables and plot type. Frames are
computed from wavefunctions, not prerecorded movies. Cloud Node numerical
parity is distinct from browser rendering QA. Browser rendering must not be
claimed from syntax checks or accessed by bypassing an existing URL policy.

New systems can reuse state, scheduling, evidence, comparison and checkpoint
logic. They MUST provide their own justified Hamiltonian, representation,
boundary, preparation, controls, observables, appropriate propagator and
benchmarks. The spectral split operator is not suitable for every system.
No autonomous discovery, scientific novelty or global solver qualification.
