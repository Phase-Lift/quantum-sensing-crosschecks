# Dynamical cross-check benchmark v0.3

**Result:** A drive-off, two-amplitude initial-velocity measurement separates the
specified GUP effective model from an actuator-only counterfeit and a nonlinear
spring counterpart. It cannot separate an exactly matching conventional kinetic
Hamiltonian. This is a simulated effective-model test, not a hardware experiment
or evidence of sub-Planck sensitivity.

Start with `CHECKPOINT.md` for the assumptions and derivation.

## Reproduce

```sh
python -m pip install -r requirements.txt
python -m unittest -v test_crosscheck.py
python oscillator_crosscheck.py --output .
```

The recorded run uses 5,000 repeats per scenario, a 64-level oscillator
representation, 32 million hypothetical destructive position readings per
experiment, plus 100,000 GHZ preparations. Sampling uses aggregated multinomial
and binomial draws; it does not individually iterate over every hypothetical
measurement. Use `--repetitions 100` for a smaller Monte Carlo check.

NumPy, SciPy and SymPy are required. Recorded versions are in
`requirements-recorded.txt` and `results.json`. Results may show tiny differences
across environments. Optional environment variables `OPENBLAS_NUM_THREADS=1`
and `OMP_NUM_THREADS=1` reduce thread overhead on some systems.

## Important boundaries

The commutator does not by itself imply the added free Hamiltonian. That dynamics
is explicitly assumed. Independent launch, readout and nuisance calibration
are required and have not been implemented. Ordinary free kinetic nonlinearity
remains exactly degenerate with the candidate parameter. The code's calibration
interval examples use hypothetical external bounds, not measured capabilities.
More repetitions do not resolve exact parameter nonidentifiability.

This stage is an extension of the preceding control benchmark; its recorded results remain separately identified.
