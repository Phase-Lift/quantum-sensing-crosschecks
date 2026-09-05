# Joint Visibility Benchmark v0.1

Status: Derived toy-model result and reproducible classical simulation of specified quantum probabilities. **Not a quantum hardware experiment, a quantum-gravity prediction, or a sub-Planck observation.**

## 1. Question and answer

Can a parameter be invisible to every individual observer, and even every pair, while remaining visible to a three-sensor network?

Yes. A standard GHZ phase family provides an explicit quantum-compatible example. The present contribution is a worked sensor-network benchmark with exact probability checks, finite-sample error calculations, noise, controls, and a nuisance-identifiability audit. The GHZ principle is established prior art, not a new theorem about gravity.

For these benchmarks, distinguishability is defined by the complete accessible joint-output distribution. Statistical processing cannot recover a distinction absent from that distribution. The present model concerns repeated probabilistic measurements, not a claim about an underlying microscopic state.

## 2. Assumptions and state preparation

There are n >= 2 sensors, each receiving one qubit of a freshly prepared state. Different trials are independent and identically distributed. Measurement bases share a calibrated phase reference. Shot identities are preserved, and visibility v is known and phase-independent.

Write |0n> = |0...0> and |1n> = |1...1>. Define

\[
\rho_\phi=\frac12\left[
|0n\rangle\langle0n|+|1n\rangle\langle1n|
+v e^{-i\phi}|0n\rangle\langle1n|
+v e^{i\phi}|1n\rangle\langle0n|
\right],\qquad 0\leq v\leq1.
\]

The nonzero eigenvalues are (1+v)/2 and (1-v)/2, so this is a valid density matrix. At v=1 it is the pure state

\[
|\psi_\phi\rangle=(|0n\rangle+e^{i\phi}|1n\rangle)/\sqrt2.
\]

For example, starting from phase zero, the unitary exp(-i phi Z_1/2) on the first qubit encodes this relative phase. This supplies an ordinary quantum encoding, **not** a mechanism by which sub-Planck structure generates phi.

## 3. Exact proper-subset blindness

For any nonempty proper subset S of the sensors, tracing out at least one qubit removes the off-diagonal coherence because <0|1>=0. Therefore

\[
\rho_S=\frac12\left(|0^{|S|}\rangle\langle0^{|S|}|+
|1^{|S|}\rangle\langle1^{|S|}|\right),
\]

independent of phi and v.

Thus every experiment confined to a fixed proper subset of these prepared probes has identical statistics for all phi. This includes collective measurements on that subset across independent copies. It does not include messages or conditioning data from the omitted sensors, changing the input-probe preparation, or adding a different physical interaction.

For equal-prior hypotheses differing only in phi, any classifier restricted to such a subset has exactly 50% optimal accuracy, for any number of copies. This is stronger than saying its mean signal is small.

## 4. Full-fleet measurement and probability law

Measure Pauli X on sensors 1 through n-1 and Pauli Y on sensor n. These are local measurements with outcomes x_i in {-1,+1}; no physically joint detector is required. Afterwards, ordinary communication allows the fleet to combine the readings.

The exact probability is

\[
P_\phi(x_1,\ldots,x_n)=2^{-n}\left[1+v\sin\phi\prod_i x_i\right].
\]

Summing over any missing x_i cancels the phase-dependent term. Every proper measurement marginal is therefore uniform. The full parity W=product_i x_i instead obeys

\[
P_\phi(W=+1)=\frac{1+v\sin\phi}{2},\qquad
\mathbb E_\phi W=v\sin\phi.
\]

The full record's phase likelihood depends only on W; parity is a sufficient statistic for phi for this selected measurement. Other details of the record should still be retained for diagnosing violations of the assumed noise model.

The classical Fisher information per trial is

\[
F_\phi=\frac{v^2\cos^2\phi}{1-v^2\sin^2\phi}.
\]

This is positive where v>0 and cos(phi) != 0 for the noisy interior model. Near phase zero, the local standard-error scale is

\[
\sigma_\phi\simeq\frac{1}{v\sqrt N}.
\]

The sine readout alone is not globally identifying: sin(phi)=sin(pi-phi). A second all-X setting yields v cos(phi); both quadratures recover phi modulo 2*pi, assuming v>0 and known phase-reference conventions. Those extra trials must be included in the resource budget. The benchmark's two hypotheses are on the local, unambiguous branch.

## 5. Numerical benchmark actually run

Parameters: n=3, v=0.8, H0: phi=0, H1: phi=0.02 radians. Visibility here means coherence contrast, not a fraction of correct trials.

The predicted correlation under H1 is

\[
C_1=0.8\sin(0.02)=0.015998933354666466.
\]

If K is the number of positive parities in N trials, then K is binomial with p0=0.5 or p1=(1+C1)/2. The code selects the exact equal-prior likelihood-ratio decision threshold and evaluates both error tails.

| Three-sensor trials | Local readings | Exact optimal binary accuracy | Monte Carlo accuracy | Any fixed proper subset |
|---:|---:|---:|---:|---:|
| 1,000 | 3,000 | 59.9833% | 59.9825% | 50% |
| 10,000 | 30,000 | 78.8135% | 78.8800% | 50% |
| 100,000 | 300,000 | 99.4293% | 99.5150% | 50% |

Monte Carlo used 20,000 independently simulated datasets under each hypothesis at each size, with seed 20260904. It sampled the exact binomial sufficient statistic; a separate raw-triple sampler checked an explicit aligned record. Results for 1,000,000 trials are also in results.json; rounded near-100% accuracy is not a claim of logically certain discrimination.

At N=100,000 the exact false-positive probability is 0.00575746 and the false-negative probability is 0.00565630. These are classification performance numbers, **not** a posterior probability that a physical theory is true and not a discovery-significance claim.

One raw 100,000-trial H1 record gave C_hat=0.01852 and

\[
\hat\phi=0.0231521\ \mathrm{rad},\qquad
\mathrm{SE}(\hat\phi)\simeq0.00395313\ \mathrm{rad}.
\]

This is compatible with the injected value 0.02 rad. In this run a zero-phase control gave parity -0.00646, and independently shuffling the third sensor's shot alignment gave parity 0.0. A single null sample can fluctuate; zero is an expectation, not a requirement on every dataset. These controls do not rule out every correlated-noise mechanism.

Exact matrix checks verified all six proper reductions for three sensors, normalized measurement probabilities, and the Born-rule formula. The maximum reduction difference and Born-formula discrepancy were both 0 in this run. The smallest computed density eigenvalue was -6.85e-32, consistent with roundoff around an exact zero; checks used a 1e-14 tolerance.

## 6. The attribution problem: seeing a phase is not identifying its cause

Suppose the observed phase is phi=theta+b, where theta is a proposed new effect and b is an ordinary unknown phase offset. Then the two models

- theta=0.02, b=0.03;
- theta=0, b=0.05

produce exactly the same full density matrix, not merely similar selected measurement statistics. No measurement of that state, fleet analysis, memory completion, or increased sample count separates them without additional assumptions or controls.

For the sine readout the Fisher information matrix for (theta,b) is proportional to

\[
\begin{pmatrix}1&1\\1&1\end{pmatrix},
\]

which has rank one. This gives a precise attribution failure even though total phase is measurable.

### Conditional control repair, also simulated

Assume two calibrated settings are available such that

\[
\phi_+=b+\theta,\qquad \phi_-=b-\theta.
\]

The setting reversal is assumed to flip theta while leaving the same b and v unchanged. Merely reversing a measurement sign flips all contributions and does not provide this discrimination. A real proposal must derive a physical response difference.

On a known local arcsine branch,

\[
\hat\theta=\frac12\left[\arcsin(\hat C_+/v)-\arcsin(\hat C_-/v)\right],
\quad
\hat b=\frac12\left[\arcsin(\hat C_+/v)+\arcsin(\hat C_-/v)\right].
\]

Using 50,000 trials per setting (the same total 100,000 trials and 300,000 local readings), theta=0.02, and b=0.03, the information matrix becomes full rank. Over 20,000 simulated experiments, the mean theta estimate was 0.01999695 rad, with empirical standard deviation 0.00392829 rad. The local Cramer-Rao standard-deviation scale is 0.00395377 rad. The small difference reflects finite Monte Carlo sampling and the local/asymptotic nature of the bound; it is not a claim of beating the bound.

No Planck-sensitive interaction with this discriminating reversal has been established here. Uncontrolled bias drift between settings or other unknown parameters can destroy the separation.

## 7. Resource accounting and fleet size

One trial consumes one n-qubit preparation and n local readings. Preserving or discarding parts of the same record is an **information-ablation test**, not a proof of superiority over the best independent sensor, interferometer, or entangled-probe strategy.

A single integrated instrument with access to the same prepared state and the same measurements reproduces the parity statistic. The result establishes the value of an accessible correlation, not a special benefit of spatial distribution.

For an illustrative independent readout-flip probability q per sensor, parity contrast is multiplied by (1-2q)^n. With total encoded phase fixed and a fixed local-reading budget R, the expected parity shift over its null standard error is

\[
(1-2q)^n\sin\phi\sqrt{R/n}.
\]

For q=0.05, phi=0.02, and R=300,000, this quantity is 6.274 for n=2, 4.610 for n=3, 2.893 for n=5, and 1.208 for n=10. Thus bigger fleets are worse in this particular fixed-total-phase noise model. This is not universal: different encodings can make the accumulated phase itself depend on n.

## 8. What connects this to the Planck question?

Nothing in the model assigns a length to phi. Calling an injected phase a sub-Planck feature does not turn the calculation into quantum gravity.

The missing physical calculation has the form

\[
\text{specified microscopic model}\longrightarrow\delta H
\longrightarrow \text{probe state/channel}\longrightarrow
P(\text{fleet data}\mid\text{model, controls, noise}).
\]

It must establish the allowed coupling, scale dependence, state preparation and measurement resources, and a response distinguishable from ordinary alternatives. A measurable effective parameter need not uniquely identify microscopic structure or its length scale.

For scale illustration only, ASSUME phi=ell/L, with L=1 m. Taking ell=0.1 times the NIST 2022 CODATA Planck length gives ell=1.616255e-36 m. Under the independent-trial noise model with v=0.8, the small-phase number of trials for an expected shift of five null standard errors is

\[
N\simeq\frac{25}{v^2(\ell/L)^2}\simeq1.49534\times10^{73}.
\]

This assumed coupling is not derived from gravity. The number is neither a forecast for all quantum-gravity models nor a universal bound on length measurements. It only demonstrates why a dimensionless phase sensitivity cannot be silently reinterpreted as sub-Planck spatial resolution. An expected five-standard-error shift also does not specify a detection-power requirement.

## 9. Stochastic fleet-blindness statement

If two models give the same distribution of the complete accessible record D, any randomized analysis represented by a parameter-independent kernel K(A|D) also has the same output distribution:

\[
P_j(A)=\int K(A\mid D)P_j(dD),\qquad j=0,1.
\]

Equality of P0(D) and P1(D) implies equality after the integral. For adaptive experiments, the equality must hold for the complete transcripts of the permitted protocols, not just for one initial measurement. Choosing genuinely different physical interventions can change the available record and its distinguishability.

## 10. Reproduction and references

Run:

```sh
python -m pip install numpy scipy
python joint_visibility.py --seed 20260904 --trials 20000 --output results.json
```

The recorded run used NumPy 2.3.5 and SciPy 1.17.0. Floating-point details or random draws may differ across environments. Exact formulas and assumptions are specified independently of numerical sampling.

Established prior art and physical context:

1. Walck and Lyons, *Only n-Qubit Greenberger-Horne-Zeilinger States Are Undetermined by Their Reduced Density Matrices*, Phys. Rev. Lett. 100, 050501 (2008). https://arxiv.org/abs/0707.4428
2. Proctor, Knott, and Dunningham, *Multiparameter Estimation in Networked Quantum Sensors*, Phys. Rev. Lett. 120, 080501 (2018). https://arxiv.org/abs/1707.06252
3. *Quantum-private distributed sensing*. https://arxiv.org/abs/2410.00970
4. Donoghue, *Quantum General Relativity and Effective Field Theory*. https://arxiv.org/abs/2211.09902
5. NIST, 2022 CODATA Planck length. https://physics.nist.gov/cgi-bin/cuu/Value?plkl=

**Conclusion:** The correlation-visibility test passes. Unique microscopic attribution and physically justified sub-Planck sensitivity have not been established. The next research target is a specified microscopic interaction and a control response that separates it from ordinary phase effects.
