# Quantum World Core-01

Executable collinear neutral-H3 nuclear wavepacket world on the MIT-licensed
QuantumModelLib LSTH ground-doublet Born–Oppenheimer surface. Not a new PES,
not interpolation of the 76 unqualified static H3 points. The six functional
modules are in `world.py`. Browser implementation uses the same sampled
analytic PES and FFT kinetic symbol, not static-data fitting or movie playback.

## Reproduce

Build the Dockerfile with the pinned BASE_IMAGE in the manual workflow, then
use `host_run.py` on an admitted Linux standard runner. No dependencies are
installed in production environments. The one manual cloud run has one CPU,
one numerical thread, 1GiB memory, no extra swap/network, 15min job ceiling.
No paid runner, artifact/cache/Packages storage, or automatic rerun.

For Python: `PhysicsEngine`, `ExperimentControl.prepare/apply/step/observe/fork/
compare/export_evidence`. `apply` currently implements only bounded momentum
kicks, with work accounting. Arbitrary laser/control fields are NOT implemented.
All energies/times/distances use explicitly recorded atomic units.

The generated `world.html` is standalone: change preparation, start/pause,
step, branch, kick, compare, export. Its heatmap is internal-coordinate
probability density. Its mean-nuclear-position diagram is illustrative, not
a classical trajectory or full three-dimensional probability reconstruction.
Browser numerical-core parity is tested in cloud Node; visual rendering QA
is separately reported, never inferred from numerical tests.

Historical E01-C/C2 NO-GO, H2 baseline NOT_QUALIFIED, numerical solver INCOMPLETE,
and original discovery NOT_TESTED remain unchanged. Scene01 original attachment
bytes are unavailable; Scene02 source/report/API and actual 76 records are
read, their failures retained, not used to construct a PES.

## Sources and physical ceiling

QuantumModelLib commit73f2ef3d322b253a4ef39e08216261bd974180f6,
`SRC/QML/H3_m.f90` and `Tests/RES_old/H3_LSTH.txt`. Original MIT notice is retained.
The scalar routines/constants are extracted unchanged, only a C ABI wrapper
is added. LSTH references Siegbahn & Liu JCP68,2457(1978), Truhlar & Horowitz
JCP68,2466(1978), DOI10.1063/1.436019. This older approximate PES is not a modern
high-accuracy electronic-structure replacement. NIST CODATA2018 proton mass
ratio1836.15267343 is the chosen nuclear BO kinetic convention, rather than
the atomic-mass convention used in some scattering literature.

Closed-box norm/energy and analytic Gaussian propagation checks address the
propagator. Upstream energies address implementation fidelity, not independent
physical accuracy. Published rounded LSTH collinear barrier0.425eV is separately
compared. CAP loss is accounted, not called product flux. Regional populations
are finite-time, finite-box observations, not converged state-resolved rates.

Free standard-runner basis: https://docs.github.com/en/billing/concepts/product-billing/github-actions
No artifact/cache upload is used; actual outputs are recovered from ordinary
logs and saved as non-sensitive Git files. Repository source history remains
immutable. No automatic ongoing hosted service is claimed.
