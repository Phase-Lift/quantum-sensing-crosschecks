# Quantum sensing cross-checks

**Status: reproducible analytical and computational benchmarks. No hardware
experiment, verified quantum-gravity effect, or sub-Planck observation is claimed.**

This project asks a practical inference question: when several quantum sensors
report a correlated signal, which physical explanations can their data actually
separate? It contains three successive benchmarks, including explicit conventional
counterexamples. The most important result is not merely a nonzero signal; it is
a calculation of what that signal does and does not identify.

## Start here

| Stage | Result within the stated model | Main limit |
|---|---|---|
| [Joint visibility](benchmarks/01_joint_visibility/CHECKPOINT.md) | A GHZ relative phase is invisible to every proper subset and visible in the complete parity record. | Established correlation principle; not evidence of microscopic length resolution. |
| [Commutator controls](benchmarks/02_commutator_controls/CHECKPOINT.md) | A specified 32-pulse sequence isolates a conditional three-sensor phase. | An ordinary nonlinear actuator gives the same modeled pulse channel. |
| [Dynamical cross-check](benchmarks/03_dynamical_crosscheck/CHECKPOINT.md) | A drive-off, two-amplitude initial-velocity contrast separates the candidate model from specified actuator and spring imitations. | An exactly matching conventional kinetic term remains indistinguishable. |

For the final stage, the independently observed effective coefficients are

\[
L=\beta+\eta,\qquad D=\beta+\kappa.
\]

Here \(\beta\) is a hypothesized deformation, \(\eta\) a conventional actuator
nonlinearity, and \(\kappa\) a conventional free-kinetic nonlinearity. Without
independently justified constraints, the two observations do not identify
\(\beta\) separately. More trials do not resolve this exact parameter redundancy.
See the final-stage checkpoint, sections 4 and 7, for the derivation and conditional
calibration intervals. Numerical inputs used for those intervals are hypothetical,
not experimentally achieved calibration bounds.

## Run the checks

Use Python 3.10 or newer in a separate environment:

```sh
python -m pip install -r requirements.txt
python run_checks.py
```

This runs the six final-stage unit tests, two small sampling checks, and the
commutator-control benchmark with a reduced repetition count. New outputs go to
`outputs/quick/`; recorded result files are not overwritten. On systems with many
BLAS threads, the runner limits BLAS/OpenMP threads to one per subprocess.

For the originally specified Monte Carlo repetition counts:

```sh
python run_checks.py --full
```

The full mode writes `outputs/full/`. A smaller run tests execution and model
consistency; it is not a replacement for the recorded high-repetition results.
Library versions and frozen parameters accompany the original numerical records.
Floating-point and random-number implementation differences can change outputs.

To check the release files for accidental changes:

```sh
python verify_files.py
```

## Scientific boundaries

The generalized-uncertainty-principle commutator is a phenomenological test
hypothesis. Its free Hamiltonian is an additional assumption, not a result implied
by the commutator alone. The polynomial simulator is not a complete microscopic
or all-orders minimal-length theory. The toy coefficients deliberately make the
responses numerically measurable; they do not demonstrate Planck-scale reach.

State preparation, control Hamiltonians, readout linearity, nuisance bounds,
dissipation and systematic errors need physical validation. Finite-basis checks
verify selected numerical states, not an exact finite-dimensional canonical algebra.

The GHZ principle, group symmetrization, commutator-loop proposals, and free
oscillator tests have established prior art. Original research references and
attribution are retained in the stage checkpoints. Novelty of these particular
benchmark compositions has not been established by an exhaustive prior-art search.

## Review and provenance

The derivations, code and documentation were developed with AI assistance. The
maintainer is responsible for checking their claims. This is not independent peer
review. See [PROVENANCE.md](PROVENANCE.md) and the packaging verification record.

Useful reviews include reproducing numerical results, finding algebraic mistakes,
constructing further conventional counterexamples, and specifying independently
calibrated experimental controls. A failed attribution test is a result to report,
not a reason to hide a counterexample.

No new distribution license was selected during packaging. See
[LICENSE_STATUS.md](LICENSE_STATUS.md) before treating this as an open-source release.
