#!/usr/bin/env python3
"""Reproducible GHZ sensor-network visibility benchmark.

This is a classical numerical simulation of specified quantum probabilities,
NOT a quantum hardware run or a test of Planck/sub-Planck physics.

Dependencies: Python >= 3.10, numpy, scipy.
Run: python joint_visibility.py --output results.json
"""
from __future__ import annotations

import argparse
from itertools import combinations, product
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.stats import binom

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int8]


def ghz_density(n: int, phase: float, visibility: float) -> ComplexArray:
    """GHZ coherence reduced by a known phase-independent visibility."""
    if not 2 <= n <= 10:
        raise ValueError("Dense verification supports 2 <= n <= 10.")
    if not np.isfinite(phase) or not 0.0 <= visibility <= 1.0:
        raise ValueError("Require finite phase and 0 <= visibility <= 1.")
    rho = np.zeros((2**n, 2**n), dtype=np.complex128)
    rho[0, 0] = rho[-1, -1] = 0.5
    rho[0, -1] = visibility * np.exp(-1j * phase) / 2
    rho[-1, 0] = np.conjugate(rho[0, -1])
    return rho


def partial_trace(rho: ComplexArray, keep: Sequence[int], n: int) -> ComplexArray:
    kept = tuple(sorted(keep))
    if len(set(kept)) != len(kept) or any(i < 0 or i >= n for i in kept):
        raise ValueError("Invalid subsystem selection.")
    traced = tuple(i for i in range(n) if i not in kept)
    order = kept + traced + tuple(i + n for i in kept) + tuple(i + n for i in traced)
    dk, dt = 2**len(kept), 2**len(traced)
    tensor = rho.reshape((2,) * (2 * n)).transpose(order).reshape(dk, dt, dk, dt)
    return np.einsum("atbt->ab", tensor)


def kron_all(matrices: Sequence[ComplexArray]) -> ComplexArray:
    out = np.array([[1.0]], dtype=np.complex128)
    for matrix in matrices:
        out = np.kron(out, matrix)
    return out


def measurement_probabilities(n: int, phase: float, visibility: float) -> tuple[FloatArray, FloatArray]:
    """Measure X on n-1 qubits and Y on the final qubit."""
    x = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    identity = np.eye(2, dtype=np.complex128)
    operators = [x] * (n - 1) + [y]
    rho = ghz_density(n, phase, visibility)
    born, analytic = [], []
    for outcomes in product((-1, 1), repeat=n):
        projector = kron_all([(identity + a * op) / 2 for a, op in zip(outcomes, operators)])
        born.append(float(np.trace(rho @ projector).real))
        analytic.append((1 + visibility * np.sin(phase) * np.prod(outcomes)) / 2**n)
    return np.array(born), np.array(analytic)


def sample_records(n: int, shots: int, phase: float, visibility: float, rng: np.random.Generator) -> IntArray:
    """Exact classical sampler for the selected measurement distribution.

    A sample is an aligned row of n local +/-1 readings. This classical
    sampler does not demonstrate entanglement or a physical coupling.
    """
    if n < 2 or shots < 1 or not 0 <= visibility <= 1:
        raise ValueError("Invalid sample parameters.")
    parity = np.where(rng.random(shots) < (1 + visibility * np.sin(phase)) / 2, 1, -1)
    first = rng.choice(np.array([-1, 1], dtype=np.int8), size=(shots, n - 1))
    last = parity * np.prod(first, axis=1)
    return np.column_stack((first, last)).astype(np.int8)


def classifier_summary(shots: int, p0: float, p1: float, trials: int, rng: np.random.Generator) -> dict:
    """Equal-prior optimal likelihood test for two specified binomial models."""
    if not 0 < p0 < p1 < 1:
        raise ValueError("Require 0 < p0 < p1 < 1.")
    log_failure_ratio = np.log1p(-p1) - np.log1p(-p0)
    log_success_ratio = np.log(p1 / p0)
    boundary = -shots * log_failure_ratio / (log_success_ratio - log_failure_ratio)
    threshold = int(np.ceil(boundary))
    false_positive = float(binom.sf(threshold - 1, shots, p0))
    false_negative = float(binom.cdf(threshold - 1, shots, p1))
    error = (false_positive + false_negative) / 2
    k0 = rng.binomial(shots, p0, size=trials)
    k1 = rng.binomial(shots, p1, size=trials)
    mc_accuracy = float(((k0 < threshold).mean() + (k1 >= threshold).mean()) / 2)
    mc_se = float(np.sqrt(((1 - false_positive) * false_positive + (1 - false_negative) * false_negative) / (4 * trials)))
    return {
        "shots": shots,
        "sensor_readings": 3 * shots,
        "decision_H1_when_positive_parity_count_at_least": threshold,
        "false_positive_probability": false_positive,
        "false_negative_probability": false_negative,
        "exact_equal_prior_accuracy": 1 - error,
        "monte_carlo_equal_prior_accuracy": mc_accuracy,
        "monte_carlo_standard_error": mc_se,
        "monte_carlo_trials_per_hypothesis": trials,
        "any_fixed_proper_subset_exact_accuracy": 0.5,
    }


def fisher_matrix(settings: Sequence[float], shots: Sequence[int], theta: float, bias: float, visibility: float) -> FloatArray:
    """Classical FI for phase_j = setting_j * theta + bias.

    Parameter order: theta, bias. Visibility and response settings are
    assumed calibrated; each setting has the stated number of fresh shots.
    """
    matrix = np.zeros((2, 2))
    for setting, count in zip(settings, shots):
        phase = setting * theta + bias
        corr = visibility * np.sin(phase)
        weight = count * (visibility * np.cos(phase))**2 / (1 - corr**2)
        derivative = np.array([setting, 1.0])
        matrix += weight * np.outer(derivative, derivative)
    return matrix


def run(seed: int, trials: int) -> dict:
    if trials < 100:
        raise ValueError("Use at least 100 Monte Carlo trials per hypothesis.")
    rng = np.random.default_rng(seed)
    n, visibility, phase0, phase1 = 3, 0.8, 0.0, 0.02
    rho0 = ghz_density(n, phase0, visibility)
    rho1 = ghz_density(n, phase1, visibility)
    marginal_differences = {}
    for k in range(1, n):
        for kept in combinations(range(n), k):
            diff = np.max(np.abs(partial_trace(rho0, kept, n) - partial_trace(rho1, kept, n)))
            marginal_differences[",".join(map(str, kept))] = float(diff)
    born_errors = []
    min_eigenvalue = 1.0
    density_trace_error = 0.0
    for phase in (phase0, phase1, -0.7, 1.2):
        rho = ghz_density(n, phase, visibility)
        min_eigenvalue = min(min_eigenvalue, float(np.linalg.eigvalsh(rho).min()))
        density_trace_error = max(density_trace_error, float(abs(np.trace(rho) - 1)))
        born, analytic = measurement_probabilities(n, phase, visibility)
        born_errors.append(float(np.max(np.abs(born - analytic))))
        assert np.all(born >= -1e-14) and np.isclose(np.sum(born), 1)
    assert max(marginal_differences.values()) < 1e-14
    assert max(born_errors) < 1e-14
    assert min_eigenvalue >= -1e-14
    assert density_trace_error < 1e-14

    corr = float(visibility * np.sin(phase1))
    p0, p1 = 0.5, (1 + corr) / 2
    classifications = [classifier_summary(shots, p0, p1, trials, rng) for shots in (1000, 10000, 100000, 1000000)]

    # Full raw-data demo and two diagnostic controls.
    shots = 100000
    records = sample_records(n, shots, phase1, visibility, rng)
    null_records = sample_records(n, shots, phase0, visibility, rng)
    shuffled = records.copy()
    shuffled[:, -1] = shuffled[rng.permutation(shots), -1]
    sample_c = float(np.prod(records, axis=1).mean())
    sample_phase = float(np.arcsin(np.clip(sample_c / visibility, -1, 1)))
    example = {
        "shots": shots,
        "single_means": records.mean(axis=0).tolist(),
        "pair_means": {f"{i},{j}": float((records[:, i] * records[:, j]).mean()) for i, j in combinations(range(n), 2)},
        "triple_correlation": sample_c,
        "expected_triple_correlation": corr,
        "estimated_phase_rad": sample_phase,
        "true_phase_rad": phase1,
        "phase_standard_error_delta_method": float(np.sqrt((1 - corr**2) / shots) / (visibility * np.cos(phase1))),
        "zero_phase_control_correlation": float(np.prod(null_records, axis=1).mean()),
        "shuffled_alignment_control_correlation": float(np.prod(shuffled, axis=1).mean()),
        "controls_note": "These detect selected pipeline errors, not all ordinary correlated noise or quantum-gravity alternatives.",
    }

    # Nuisance non-identifiability, then an explicitly assumed response reversal.
    theta, bias, total_shots = 0.02, 0.03, 100000
    single_fi = fisher_matrix([1], [total_shots], theta, bias, visibility)
    reversed_fi = fisher_matrix([1, -1], [total_shots // 2] * 2, theta, bias, visibility)
    candidate_a = ghz_density(3, theta + bias, visibility)
    candidate_b = ghz_density(3, 0.0 + (theta + bias), visibility)
    assert np.max(np.abs(candidate_a - candidate_b)) == 0.0
    assert np.linalg.matrix_rank(single_fi) == 1
    assert np.linalg.matrix_rank(reversed_fi) == 2
    n_each = total_shots // 2
    c_plus = visibility * np.sin(bias + theta)
    c_minus = visibility * np.sin(bias - theta)
    count_plus = rng.binomial(n_each, (1 + c_plus) / 2, trials)
    count_minus = rng.binomial(n_each, (1 + c_minus) / 2, trials)
    estimate_plus = np.arcsin(np.clip((2 * count_plus / n_each - 1) / visibility, -1, 1))
    estimate_minus = np.arcsin(np.clip((2 * count_minus / n_each - 1) / visibility, -1, 1))
    theta_estimates = (estimate_plus - estimate_minus) / 2
    bias_estimates = (estimate_plus + estimate_minus) / 2
    nuisance = {
        "assumption": "Known setting reversal flips theta but leaves the same additive bias b and visibility unchanged; all phases stay on the local arcsin branch.",
        "theta_rad": theta,
        "bias_rad": bias,
        "total_shots": total_shots,
        "total_sensor_readings": 3 * total_shots,
        "single_setting_fisher_matrix": single_fi.tolist(),
        "single_setting_rank": int(np.linalg.matrix_rank(single_fi)),
        "exact_degenerate_single_setting_models": [{"theta": theta, "bias": bias}, {"theta": 0.0, "bias": theta + bias}],
        "two_setting_fisher_matrix": reversed_fi.tolist(),
        "two_setting_rank": int(np.linalg.matrix_rank(reversed_fi)),
        "local_CR_standard_deviation_theta": float(np.sqrt(np.linalg.inv(reversed_fi)[0, 0])),
        "monte_carlo_theta_mean": float(theta_estimates.mean()),
        "monte_carlo_theta_sd": float(theta_estimates.std(ddof=1)),
        "monte_carlo_bias_mean": float(bias_estimates.mean()),
        "warning": "No physical mechanism providing this discriminating reversal has been derived for Planck-scale physics.",
    }

    # Illustrative fleet-size audit. Keep TOTAL encoded phase and total read
    # budget fixed. This is not an optimized-metrology or universal bound.
    size_audit = []
    total_readings, q = 300000, 0.05
    for fleet_size in (2, 3, 5, 10):
        used_shots = total_readings // fleet_size
        v = (1 - 2 * q)**fleet_size
        c = v * np.sin(phase1)
        size_audit.append({
            "fleet_size": fleet_size,
            "shots": used_shots,
            "independent_readout_error_per_sensor": q,
            "effective_visibility": v,
            "mean_parity_shift_over_null_standard_error": float(c * np.sqrt(used_shots)),
        })

    # This weak-coupling example is ASSUMED, not derived from gravity.
    planck_length = 1.616255e-35  # m, NIST 2022 CODATA central value
    hypothetical_length = 0.1 * planck_length
    illustrative_required_shots = (5 / (visibility * hypothetical_length / 1.0))**2
    return {
        "status": "TOY_MODEL_ONLY_NOT_A_PLANCK_SCALE_OBSERVATION",
        "seed": seed,
        "model": {
            "n": n, "visibility": visibility,
            "H0_phase_rad": phase0, "H1_phase_rad": phase1,
            "measurement": "X tensor X tensor Y",
            "probability": "P(x1,x2,x3|phi) = (1+v*sin(phi)*x1*x2*x3)/8",
            "assumptions": ["Independent identical state preparations", "Calibrated phase reference and visibility", "Aligned shot identities", "No adversarial or unmodeled correlated noise", "The phase is an input, not a derived sub-Planck signal"],
        },
        "exact_checks": {
            "proper_reduced_state_max_differences": marginal_differences,
            "max_Born_rule_vs_closed_form_error": max(born_errors),
            "minimum_density_eigenvalue": min_eigenvalue,
            "maximum_density_trace_error": density_trace_error,
            "checks_passed": True,
        },
        "classification": classifications,
        "one_raw_record_demo": example,
        "nuisance_identifiability": nuisance,
        "fleet_size_fixed_read_budget_audit": size_audit,
        "hypothetical_linear_length_coupling": {
            "assumed_relation": "phi = ell / L; NOT a gravity prediction",
            "L_m": 1.0,
            "ell_m": hypothetical_length,
            "ell_over_planck_length": 0.1,
            "visibility": visibility,
            "N_for_expected_shift_5_null_standard_errors": illustrative_required_shots,
            "note": "Approximate small-phase, iid shot-noise scaling. An expected 5-sigma shift is not a specified detection power, and this is not a universal length bound.",
        },
        "resource_limit": "A single integrated instrument with the same state and measurements reproduces the joint statistic. Proper-subset tests are ablations, not a comparison with the best alternative sensor design.",
        "provenance": {
            "numpy_version": np.__version__,
            "references": [
                "https://arxiv.org/abs/0707.4428",
                "https://arxiv.org/abs/1707.06252",
                "https://arxiv.org/abs/2410.00970",
                "https://arxiv.org/abs/2211.09902",
                "https://physics.nist.gov/cgi-bin/cuu/Value?plkl",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--trials", type=int, default=20000)
    parser.add_argument("--output", type=Path, default=Path("results.json"))
    args = parser.parse_args()
    result = run(args.seed, args.trials)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("TOY MODEL ONLY: no Planck-scale observation or quantum hardware run.")
    print("Exact state/marginal/Born-rule checks:", result["exact_checks"]["checks_passed"])
    for row in result["classification"]:
        print(f"N={row['shots']:>7}: exact={row['exact_equal_prior_accuracy']:.8%}, MC={row['monte_carlo_equal_prior_accuracy']:.8%}, proper subset=50%")
    print("Raw-record example:", json.dumps(result["one_raw_record_demo"], indent=2))
    print("Nuisance check:", json.dumps(result["nuisance_identifiability"], indent=2))
    print("Results:", args.output.resolve())


if __name__ == "__main__":
    main()
