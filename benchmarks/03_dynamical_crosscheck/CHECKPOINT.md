# Dynamical cross-check checkpoint v0.3

**Status:** Analytical effective-model derivation, symbolic checks, finite-basis quantum evolution, and sampling from calculated measurement probabilities. **Not a hardware experiment, a derivation from an underlying geometric theory, a full quantum-gravity theory, or a sub-Planck observation.**

## Executive result

The next step after v0.2 is not to add observers or repeat the same pulse experiment. It is to add an independently calibrated, drive-off measurement whose parameter dependence differs from the nonlinear actuator's response.

This checkpoint supplies that measurement: a two-amplitude, opposite-launch **initial-position-velocity contrast**. Under the specified free-evolution completion, it distinguishes the proposed deformation from both an actuator-only counterfeit and an actuator-plus-nonlinear-spring counterfeit. Measuring oscillator energy levels alone would not distinguish the second counterfeit: the relevant two polynomial Hamiltonians are exactly isospectral.

A conventional momentum-dependent free Hamiltonian remains an exact counterfeit. The resulting experiment identifies two effective coefficients, not a uniquely attributable microscopic cause:

\[
L=\beta+\eta,\qquad D=\beta+\kappa.
\]

Here beta is the candidate deformation, eta is an actuator nonlinearity, and kappa is a conventional kinetic nonlinearity. The ordinary quartic spring coefficient lambda does not enter the ideal initial-velocity contrast.

**Decision:** the cross-check passes as conditional mechanism discrimination. Exclusive microscopic attribution and physical Planck-scale sensitivity remain unestablished. The deliverable includes a calibration-aware interval routine that will not exclude beta=0 simply by assuming unknown nuisances away.

## 1. Provenance and the hypothesis added in this step

The supplied `../02_commutator_controls/CHECKPOINT.md` derived a 32-pulse three-qubit loop. At its benchmark pulse strengths, the GHZ phase is

\[
\Phi=8(\beta+\eta).
\]

That checkpoint explicitly did not specify oscillator free dynamics and did not derive its phenomenological commutator from an independently specified underlying geometric theory. We preserve those limitations.

For a **conditional dynamical completion**, assume

\[
[X,P]=i(1+\beta P^2),\qquad
H_{\rm candidate}/(\hbar\omega)=(P^2+X^2)/2.
\]

Using the local canonical representation

\[
X=x,\quad P=p+\beta p^3/3+O(\beta^2),\quad [x,p]=i,
\]

gives

\[
h_{\rm candidate}=\frac{x^2+p^2}{2}+\frac{\beta}{3}p^4+O(\beta^2).
\]

The commutator does **not** imply this Hamiltonian. It is an additional choice, comparable to a familiar phenomenological completion in existing oscillator tests [R1]. The relationship between minimal-length kinematics, free dynamics and relativity requirements is nontrivial [R2]. None of the results below apply automatically to every GUP model.

The actual simulator uses the ordinary canonical **polynomial effective Hamiltonian**

\[
\boxed{h=\frac{x^2+p^2}{2}+\frac{\beta+\kappa}{3}p^4+
\frac{\lambda}{3}x^4.}
\]

It is exact within that stipulated effective model. Its interpretation as the selected GUP is perturbative. The tested quartic coefficients are nonnegative, so these particular polynomial examples are stable. Time is dimensionless, tau=omega*t; all coefficient and quadrature values below are dimensionless.

The pulse coupling is retained as

\[
P_{\rm pulse}=p+\frac{\beta+\eta}{3}p^3.
\]

The actuator coefficient eta is present only while the suspect pulse drive is on. The kinetic coefficient kappa remains during free evolution. This is a classification of the specified model terms, not a claim that an actual instrument necessarily contains those exact nonlinearities.

## 2. Why energy spectroscopy alone is not the right next measurement

Consider

\[
h_p=h_0+\frac{a}{3}p^4,\qquad
h_x=h_0+\frac{a}{3}x^4,\qquad h_0=(x^2+p^2)/2.
\]

A quarter-turn oscillator phase-space rotation exchanges x and p up to sign and preserves h0. Consequently it maps hp into hx by a unitary conjugation. The two Hamiltonians have the **same complete energy spectrum**. This is exact for the polynomial models, not just a first-order perturbation result.

Their position dynamics with the same physically specified preparation and readout need not coincide. A unitary equivalence of energy operators does not license rotating the physical meaning of the measured position or the prepared state without accounting for it.

For an ordinary harmonic number state,

\[
\langle n|x^4|n\rangle=\langle n|p^4|n\rangle
=\frac34(2n^2+2n+1).
\]

First-order perturbation theory therefore gives

\[
E_n/(\hbar\omega)=n+\tfrac12+
\frac{\beta+\kappa+\lambda}{4}(2n^2+2n+1)+O(c^2),
\]

where c denotes the small nonlinear coefficients. Adjacent-transition spacing differences measure the same summed coefficient at this order.

Thus the counterfeit beta=0, eta=a, lambda=a, kappa=0 matches both the ideal pulse-loop signal and the exact energy spectrum of beta=a, eta=kappa=lambda=0. A pulse-plus-spectrum comparison alone would wrongly appear to corroborate the candidate.

Existing experimental work also reports structural oscillator nonlinearity as a limitation to attributing nonlinear dynamics to proposed gravity effects [R1]. We do not claim novelty for oscillator nonlinearity as a diagnostic or for the general spectroscopic ambiguity.

## 3. The drive-off initial-velocity contrast

Write A=beta+kappa. During drive-off evolution, the Heisenberg equation for the measured position gives

\[
\frac{dx}{d\tau}=i[h,x]=p+\frac{4A}{3}p^3.
\]

The x-only potential commutes with x. It therefore does not appear in the instantaneous velocity operator. It will affect the subsequent trajectory through momentum evolution.

Prepare identical base packets with zero mean canonical p, equal variance V, and matched preparation errors. Launch them with positive and negative mean canonical momentum q using a separately calibrated **x-coupled** operation, not the suspect nonlinear p actuator. For a Gaussian packet,

\[
v(q)=\left.\frac{d\langle x\rangle_q}{d\tau}\right|_0
=q+\frac{4A}{3}(q^3+3qV).
\]

More generally, an identical third central moment contributes an extra constant to v(+q) and v(-q) that cancels in the paired difference. Even an unknown constant coefficient c2 in the quadratic kinetic term can be allowed:

\[
R(q)=\frac{v(+q)-v(-q)}{2q}
=c_2+\frac{4A}{3}(q^2+3V).
\]

Use two different nonzero launch amplitudes q1 and q2:

\[
\boxed{D=\frac{3\,[R(q_2)-R(q_1)]}{4(q_2^2-q_1^2)}=A=\beta+\kappa.}
\]

This contrast cancels the common quadratic coefficient and the common packet variance. It does not require neglecting the wavepacket's quantum spread. It also eliminates the ideal quartic spring's direct contribution to this initial slope.

For q1=1, q2=2, the contrast simplifies to D=(R(2)-R(1))/4. For the pure candidate beta=0.0025 and V=1/2,

\[
v(1)=1.0083333333,\qquad v(2)=2.0366666667,\qquad D=0.0025.
\]

For the actuator-only or actuator-plus-quartic-spring counterpart, D=0 in the exact initial-slope limit.

### Why this is independent only under explicit calibration conditions

The variable p is the canonical variable in the stipulated representation, not an automatic identification with physical P. The preparation model assumes a known x-coupled phase gradient that translates p. Realizing and calibrating that operation independently is an uncompleted hardware requirement.

The example mu(q)=q+zeta*q^3 for the actual launch momentum provides a concrete warning. In an otherwise ordinary harmonic oscillator, incorrectly treating q as the true momentum gives D_apparent=3*zeta/4. Choosing zeta=0.0033333333 fakes D=0.0025 in this slope contrast. The fractional launch errors would be 0.3333% and 1.3333% at q=1 and q=2. This counterfeit is for the initial-slope statistic, not for the entire position trajectory.

Equal packet widths, launch amplitudes, position-readout linearity and timebase calibration must be checked independently. Constant additive readout offsets can cancel in mirrored launches; arbitrary readout nonlinearity or drift does not.

## 4. Conditional recovery and the remaining indistinguishability

If independent physical information establishes kappa=0, then

\[
\beta=D,\qquad\eta=L-D.
\]

This resolves the v0.2 actuator degeneracy **within that restricted model**.

Without that constraint, the transformation

\[
\beta\mapsto\beta+d,\quad
\eta\mapsto\eta-d,\quad
\kappa\mapsto\kappa-d
\]

leaves the pulse generator, free Hamiltonian and chosen measurement/preparation operators unchanged. Every allowed adaptive protocol built from those ingredients then has the same outcome distribution. This is stronger than merely matching one fitted phase.

Including the first-order spectroscopic coefficient S gives the response matrix

\[
\begin{pmatrix}L\\D\\S\end{pmatrix}=
\begin{pmatrix}1&1&0&0\\1&0&1&0\\1&0&1&1\end{pmatrix}
\begin{pmatrix}\beta\\\eta\\\kappa\\\lambda\end{pmatrix}.
\]

Its rank is three for four unknowns, with null direction (1,-1,-1,0). Fixing kappa externally produces full rank for beta, eta and lambda at the retained order. More precision in the same unconstrained channels does not remove the null direction.

The explicit complete conventional counterfeit is beta=0, eta=0.0025, kappa=0.0025, lambda=0. Its pulse and free generators equal those of beta=0.0025, eta=kappa=lambda=0 exactly in the simulated polynomial model.

This does not prove every microscopic theory is indistinguishable from every conventional model. It says this particular low-energy parametrization cannot uniquely identify its microscopic origin without additional physically justified constraints. No claim is made that a real apparatus has the specified kappa term at that magnitude.

## 5. Numerical experiment actually run

The program performs unitary state-vector evolution, rather than feeding an assumed slope into a noise generator. It diagonalizes the free Hamiltonian in a truncated Fock basis, computes the evolved position-measurement probabilities, and samples those probabilities with multinomial counts. Fleet parities are sampled with the previous binomial readout probability. These are simulated quantum probabilities, not quantum hardware data.

Per simulated experiment:

- Four launch settings: q=+1,-1,+2,-2, each from the same reference Gaussian with V=1/2.
- Eight positive evolution times tau=0.075,0.15,...,0.6, each on fresh preparations.
- One million destructive position readings per launch/time setting: **32 million position preparations/readings**.
- An additional 100,000 GHZ preparations split equally between bias signs: 300,000 local qubit readings. Readout visibility is 0.8 and a shared phase offset of 0.03 radians is included.

No continuous simultaneous monitoring of a single quantum oscillator is assumed. Independent fresh preparations avoid silently ignoring measurement backaction. Hardware calibration, state-preparation time, losses and dissipative overhead are not accounted for. No fixed-resource precision advantage over established sensing strategies is claimed.

The initial slope is estimated from an odd polynomial t,t^3,t^5 fitted to the positive-time position trajectory. The odd symmetry holds for the selected even Hamiltonian and chosen preparations. It is not a general-purpose fit in the presence of damping, asymmetric states or time-dependent forces.

The calculation uses 5,000 repeated simulated experiments for each of six scenarios. The reported spreads are **standard deviations of individual experiments**, not standard errors of the Monte Carlo means or probabilities that a theory is true.

| Scenario | True beta | True eta | True kappa | True lambda | Mean recovered D | SD of D |
|---|---:|---:|---:|---:|---:|---:|
| Null | 0 | 0 | 0 | 0 | 0.000002737 | 0.000590738 |
| Specified deformation | 0.0025 | 0 | 0 | 0 | 0.002489156 | 0.000593485 |
| Actuator only | 0 | 0.0025 | 0 | 0 | 0.000000641 | 0.000587619 |
| Actuator plus nonlinear spring | 0 | 0.0025 | 0 | 0.0025 | -0.000002768 | 0.000582893 |
| Complete conventional counterfeit | 0 | 0.0025 | 0.0025 | 0 | 0.002518003 | 0.000591431 |
| Mixed model with nonlinear spring | 0.0015 | 0.0010 | 0 | 0.01 | 0.001504938 | 0.000586931 |

The specified deformation and complete counterfeit have exactly the same theoretical distributions. Their slightly different sample means are Monte Carlo fluctuations, not a way to distinguish them.

For the mixed model, applying the **valid in that scenario** condition kappa=0 recovered

\[
\widehat\beta=0.00150494\quad(\mathrm{SD}=0.00058693),
\]

\[
\widehat\eta=0.00099207\quad(\mathrm{SD}=0.00076449).
\]

The inputs were beta=0.0015 and eta=0.0010. The lambda=0.01 ordinary spring term did not create a leading false kinetic signal in the contrast.

## 6. Independent numerical and approximation checks

**Dimension convergence.** For the q=2 pure-candidate trajectory over tau=0..2, dimensions 32,48,64 and 96 agreed in mean position to within 1.8e-14 relative to the 96-dimensional result. The independently evaluated commutator slope agreed with its analytical value to 4.5e-16. These are numerical consistency checks on selected low-energy states, not experimental uncertainties or a claim that finite matrices satisfy the exact canonical algebra everywhere.

**Spectrum trap.** For coefficient 0.0025, the lowest twelve eigenvalues of the p^4 and x^4 Hamiltonians agreed to 3.2e-14. The quarter-rotation operator relation was verified separately. The exact lowest transition-spacing difference was 0.0024738831, not the first-order value 0.0025. Higher-order corrections matter to high-precision spectroscopy even in the polynomial model.

**Finite-time bias.** At the selected tau_max=0.6, pure-deformation slope extrapolation had an infinite-shot coefficient bias -4.8243e-7. The actuator-plus-spring case had bias +1.3979e-7; the candidate-plus-lambda=0.01 check had bias +2.2861e-7. These are much smaller than this benchmark's statistical errors. Smaller time windows reduced the biases but increased statistical uncertainty. The check does not bound every unknown apparatus effect.

**Second-order interpretation.** The next local term in the candidate tangent representation yields

\[
h=h_0+\frac\beta3p^4+\frac{17\beta^2}{90}p^6+\cdots.
\]

For Gaussian packets, the two-amplitude estimator becomes

\[
D=\beta+\frac{17\beta^2}{20}(q_1^2+q_2^2+10V)+\cdots.
\]

At this toy parameter it is 0.002553125, a **2.125% correction**. This term was checked with matrices and independently by moment algebra. The leading inference is not an all-orders GUP result, and no global tangent-operator domain construction is attempted.

**Six executable tests** cover the algebra/rank, the complete-counterfeit generator equality, spring versus kinetic initial slopes, time-window extrapolation, absence of attribution with unbounded nuisances, and conditional interval behavior. All passed in the recorded run.

## 7. Calibration-aware inference rather than assumed attribution

If an independent estimate of kappa is available,

\[
\widehat\beta=\widehat D-\widehat\kappa.
\]

For independent unbiased errors in this same coefficient convention,

\[
\sigma_\beta^2=\sigma_D^2+\sigma_\kappa^2.
\]

Thus improving position statistics alone cannot reduce the error below the independent kappa-calibration uncertainty. Correlated errors require the covariance term; unknown bias requires a separate bound, not this variance formula.

The bundle also includes a conservative interval routine. Given independently justified hard bounds |eta|<=e and |kappa|<=k, two asymptotically normal channel errors, and a separately justified finite-window/model bias budget, it forms a simultaneous statistical rectangle and intersects the resulting allowed beta intervals. It uses z=2.2414027 for a Bonferroni 95% two-channel construction, not two separate 95% intervals mislabeled as jointly 95%.

Illustration using the **expected toy readouts** L=D=0.0025, this run's analytical standard errors, and an **assumed** additional drift-bias bound 1e-6:

| External nuisance assumptions | Allowed beta interval | Can zero be excluded within those assumptions? |
|---|---|---|
| No independent bounds | Unbounded | No |
| Both absolute bounds 0.003 | [-0.00160775, 0.00660775] | No |
| Both absolute bounds 0.0001 | [0.00129225, 0.00370775] | Yes, conditionally |

**These calibration bounds are invented test inputs, not measured apparatus capabilities.** The interval program does not establish them. The confidence statement is conditional on its statistical approximations and all error/bias bounds actually holding. An empty intersection signals inconsistency between data, model and bounds; the code does not call it a discovery. Restricting beta>=0 is an optional additional GUP hypothesis, not imposed in these examples.

A nonzero beta after justified calibration would reject beta=0 inside this chosen model class. It would still not uniquely reconstruct sub-Planck structure or identify a unique ultraviolet theory.

## 8. What remains physically unestablished

No real oscillator, independent x-coupled launch reference, verified p-coupled actuator, or position readout has been engineered or measured here. No experimentally justified kappa bound is supplied. No finite-duration pulse, damping, loss, arbitrary correlated noise, drift, or composite-body microscopic derivation is completed.

The original toy coefficient beta=0.0025 is deliberately large. The prior checkpoint's illustrative mass/frequency conversion places a Planck-order beta0=1 at beta approximately 1.5564e-41, not 0.0025. The new model's standard error near 5.9e-4 is therefore nowhere near a demonstrated Planck-scale measurement. The previous physical conversion itself assumes a particular composite-oscillator interpretation. No new length resolution is established in v0.3.

The specified dynamics are a phenomenological test branch, **not a consequence derived from an independently specified underlying geometric theory**. Free oscillator evolution as a probe is established prior art [R1], and ordinary higher-order interactions are an established issue in quantum-optical proposals [R3]. Novelty of this particular cross-check construction has not been established by an exhaustive prior-art search.

## 9. Operational decision

For this model, freeze the claim at **independently cross-checked effective nonlinearity** unless calibration or independently derived physical laws constrain the remaining conventional kinetic term. The relevant advance is that the benchmark now distinguishes some false causes and explicitly refuses to distinguish an exactly equivalent cause.

A stronger physics claim requires a specified apparatus and physically justified nuisance bounds, or a microscopic model that predicts an additional distinguishable response. More observers or more samples of the same generators cannot resolve the demonstrated parameter redundancy.

## 10. Reproduction and contents

Install dependencies and run:

```sh
python -m pip install -r requirements.txt
python -m unittest -v test_crosscheck.py
python oscillator_crosscheck.py --output .
```

`--repetitions 100` is useful for a smaller sampling check; it is not the recorded 5,000-repetition run. On systems with many BLAS threads, limiting `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` can avoid unnecessary overhead. The seed is the arbitrary fixed integer 20260906; it is not an assertion about the date of the run.

Files include the full source, six tests, frozen parameters, package versions, JSON results, trajectory and recovery CSVs, this derivation, test output and SHA-256 hashes. CSV files contain plain data, not an Excel workbook. No hardware measurements or research-independent replication are represented.

## References

[R1] M. Bawaj et al., *Probing deformed commutators with macroscopic harmonic oscillators*, Nature Communications 6, 7503 (2015). https://www.nature.com/articles/ncomms8503  
Prior experimental free-evolution strategy, additional Hamiltonian assumptions, and explicitly acknowledged structural-nonlinearity attribution limits.

[R2] P. Bosso, G. Fabiano, D. Frattulillo and F. Wagner, *Fate of Galilean relativity in minimal-length theories*, Physical Review D 109, 046016 (2024), arXiv:2307.12109. https://arxiv.org/abs/2307.12109  
Why a commutator or a minimal-length parameter alone does not specify physical dynamics. This work is not being cited as validating the phenomenological completion used here.

[R3] S. P. Kumar and M. B. Plenio, *Quantum-optical tests of Planck-scale physics*, Physical Review A 97, 063855 (2018), arXiv:1708.05659. https://arxiv.org/abs/1708.05659  
Ordinary higher-order interactions, false-positive concerns, calibration and sensitivity costs.

[R4] A. Kempf, G. Mangano and R. B. Mann, *Hilbert Space Representation of the Minimal Length Uncertainty Relation*, Physical Review D 52, 1108 (1995), arXiv:hep-th/9412167. https://arxiv.org/abs/hep-th/9412167  
Established minimal-length/GUP background; the low-order canonical polynomial simulator is not a replacement for its full operator construction.

Upstream project source: the preceding `../02_commutator_controls/CHECKPOINT.md` and its accompanying benchmark files in this repository. Our new analytic and numerical results are documented above and in `results.json` rather than attributed to these external papers.
