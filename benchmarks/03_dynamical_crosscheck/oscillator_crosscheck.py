#!/usr/bin/env python3
"""Sensor-network v0.3: drive-off dynamical cross-check and explicit identifiability audit.

A conditional effective-Hamiltonian calculation, not hardware data or evidence
of quantum gravity. H/(hbar omega)=(x^2+p^2)/2+(beta+kappa)p^4/3+lam*x^4/3.
The prior ideal pulse sequence measures L=beta+eta (phase=8L here).
The first-order GUP interpretation requires an additional free-Hamiltonian
assumption. The polynomial effective Hamiltonian is the actual simulated model.

The simulator evolves oscillator state vectors and samples finite-basis
position measurement probabilities using multinomial draws. It does NOT
insert independent Gaussian noise into analytically prescribed slopes.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.stats import norm
import sympy as sp


@dataclass(frozen=True)
class Config:
    dimension: int = 64
    repetitions: int = 5000
    position_shots_per_setting: int = 1_000_000
    fleet_trials: int = 100_000
    visibility: float = .8
    readout_offset: float = .03
    seed: int = 20260906
    tmax: float = .6
    time_points: int = 8
    momenta: tuple[float, float] = (1., 2.)


@dataclass(frozen=True)
class Case:
    name: str
    beta: float = 0.
    eta: float = 0.
    kappa: float = 0.
    lam: float = 0.

    @property
    def kinetic(self) -> float:
        return self.beta+self.kappa

    @property
    def loop(self) -> float:
        return self.beta+self.eta


CASES = [
    Case('null'),
    Case('specified_deformation', beta=.0025),
    Case('actuator_only', eta=.0025),
    Case('actuator_plus_nonlinear_spring', eta=.0025, lam=.0025),
    Case('conventional_complete_mimic', eta=.0025, kappa=.0025),
    Case('mixed_with_nonlinear_spring', beta=.0015, eta=.0010, lam=.01),
]


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def operators(dim: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    an=np.diag(np.sqrt(np.arange(1,dim)), 1)
    x=(an+an.T)/np.sqrt(2)
    p=(an-an.T)/(1j*np.sqrt(2))
    # Exact projected oscillator H0; avoids the artificial upper-edge defect
    # obtained by multiplying truncated x,p in the quadratic part.
    h0=np.diag(np.arange(dim)+.5)
    return x,p,h0


def packet(dim: int, momentum: float) -> np.ndarray:
    """Reference Gaussian with <x>=0, <p>=momentum, Var(p)=1/2.

    Independently calibrated x-coupled preparation is assumed. In particular,
    the suspect nonlinear momentum actuator is not used as its own calibrator.
    """
    a=1j*momentum/np.sqrt(2)
    psi=np.empty(dim,complex)
    psi[0]=np.exp(-abs(a)**2/2)
    for n in range(1,dim):
        psi[n]=psi[n-1]*a/np.sqrt(n)
    return psi/np.linalg.norm(psi)


def hamiltonian(case: Case, dim: int, second_order: bool=False) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    x,p,h0=operators(dim)
    h=h0+case.kinetic*np.linalg.matrix_power(p,4)/3+case.lam*np.linalg.matrix_power(x,4)/3
    if second_order:
        check(case.eta==case.kappa==case.lam==0., 'Second-order test is pure deformation only.')
        # P=tan(sqrt(beta)*p)/sqrt(beta), locally:
        # H=P^2/2+x^2/2=H0+beta p^4/3+17 beta^2 p^6/90+...
        h=h+17*case.beta**2*np.linalg.matrix_power(p,6)/90
    return h,x,p


def trajectory(case: Case, dim: int, momentum: float, times: np.ndarray,
               second_order: bool=False) -> dict[str,Any]:
    h,x,p=hamiltonian(case,dim,second_order)
    e,u=eigh(h)
    q,qvec=eigh(x)
    psi=packet(dim,momentum)
    states=u@((u.conj().T@psi)[:,None]*np.exp(-1j*e[:,None]*times))
    probs=np.abs(qvec.conj().T@states)**2
    probs/=probs.sum(axis=0,keepdims=True)
    means=q@probs
    variances=q*q@probs-means**2
    edge=np.max(np.sum(abs(states[-8:])**2,axis=0))
    initial_slope=float(np.vdot(psi,1j*(h@x-x@h)@psi).real)
    return dict(nodes=q, probabilities=probs, means=means, variances=variances,
                initial_slope=initial_slope, edge_probability=float(edge),
                eigenvalues=e, states=states)


def slope_weights(c: Config, tmax: float | None=None) -> tuple[np.ndarray,np.ndarray]:
    maximum=c.tmax if tmax is None else tmax
    times=np.linspace(maximum/c.time_points,maximum,c.time_points)
    # The selected time-reversal/parity-symmetric model with these packets has
    # an odd position trajectory. Fit t,t^3,t^5 using only nonnegative times.
    # Damping or symmetry-breaking preparation requires a different fit.
    design=np.column_stack([(times/maximum)**power for power in (1,3,5)])
    weights=np.linalg.pinv(design)[0]/maximum
    return times,weights


def drift_analysis(case: Case, c: Config, tmax: float | None=None,
                   dim: int | None=None) -> dict[str,Any]:
    times,weights=slope_weights(c,tmax)
    dimension=dim or c.dimension
    records=[]
    ratios=[]
    ratio_variances=[]
    for momentum in c.momenta:
        traces=[trajectory(case,dimension,s*momentum,times) for s in (1,-1)]
        ratio=float(weights@(traces[0]['means']-traces[1]['means'])/(2*momentum))
        ratio_var=float(np.dot(weights**2,traces[0]['variances']+traces[1]['variances'])/
                        (4*momentum**2*c.position_shots_per_setting))
        ratios.append(ratio)
        ratio_variances.append(ratio_var)
        for sign,trace in zip((1,-1),traces):
            records.append(dict(momentum=momentum,sign=sign,trace=trace))
    factor=3/(4*(c.momenta[1]**2-c.momenta[0]**2))
    fitted=factor*(ratios[1]-ratios[0])
    se=factor*np.sqrt(sum(ratio_variances))
    return dict(times=times,weights=weights,records=records,
                ratio_values=ratios,estimate=float(fitted),
                predicted_se=float(se),bias=float(fitted-case.kinetic))


def summarize(a: np.ndarray, target: float, se: float) -> dict[str,float]:
    return dict(mean=float(a.mean()),standard_deviation=float(a.std(ddof=1)),
                target=float(target),analytic_standard_error=float(se),
                coverage_95=float(np.mean(abs(a-target)<=1.959963984540054*se)))


def sample_case(case: Case,c: Config,rng: np.random.Generator) -> dict[str,Any]:
    analysis=drift_analysis(case,c)
    ratios={momentum:np.zeros(c.repetitions) for momentum in c.momenta}
    for record in analysis['records']:
        momentum,sign,trace=record['momentum'],record['sign'],record['trace']
        for j,weight in enumerate(analysis['weights']):
            # Counts of eigenvalues of the finite-basis position operator.
            counts=rng.multinomial(c.position_shots_per_setting,
                                    trace['probabilities'][:,j],size=c.repetitions)
            means=counts@trace['nodes']/c.position_shots_per_setting
            ratios[momentum]+=sign*weight*means/(2*momentum)
    factor=3/(4*(c.momenta[1]**2-c.momenta[0]**2))
    drift=factor*(ratios[c.momenta[1]]-ratios[c.momenta[0]])
    n=c.fleet_trials//2
    phase=8*case.loop
    corr=[]
    phase_vars=[]
    for sign in (1,-1):
        phi=c.readout_offset+sign*phase
        truth=c.visibility*np.sin(phi)
        k=rng.binomial(n,(1+truth)/2,size=c.repetitions)
        corr.append(2*k/n-1)
        phase_vars.append((1-truth**2)/(n*c.visibility**2*np.cos(phi)**2))
    arcs=[np.arcsin(np.clip(t/c.visibility,-1,1)) for t in corr]
    loop=(arcs[0]-arcs[1])/16
    loop_se=np.sqrt(sum(phase_vars))/16
    # Reports are conditional on kappa=0. The mimic intentionally violates it.
    eta_hat=loop-drift
    eta_se=np.sqrt(loop_se**2+analysis['predicted_se']**2)
    return dict(case=asdict(case),
                estimated_kinetic=summarize(drift,case.kinetic,analysis['predicted_se']),
                estimated_loop=summarize(loop,case.loop,float(loop_se)),
                conditional_eta=summarize(eta_hat,case.eta-case.kappa,float(eta_se)),
                estimator_bias_at_infinite_shots=analysis['bias'],
                kappa_zero_assumption_correct=(case.kappa==0),
                false_beta_detection_rate_under_kappa_zero=float(np.mean(abs(drift)>1.959963984540054*analysis['predicted_se'])))


def algebra_checks() -> dict[str,Any]:
    A,p1,p2,V,mu3,c2=sp.symbols('A p1 p2 V mu3 c2',real=True)
    def velocity(q: sp.Expr) -> sp.Expr:
        return c2*q+sp.Rational(4,3)*A*(q**3+3*q*V+mu3)
    def ratio(q: sp.Expr) -> sp.Expr:
        return (velocity(q)-velocity(-q))/(2*q)
    expr=sp.factor(sp.Rational(3,4)*(ratio(p2)-ratio(p1))/(p2**2-p1**2))
    check(expr==A,'Two-amplitude contrast did not isolate kinetic quartic.')
    # Vector order beta, eta, kappa, lambda. Spectrum row is FIRST ORDER.
    J=sp.Matrix([[1,1,0,0],[1,0,1,0],[1,0,1,1]])
    null=J.nullspace()
    gauge=sp.Matrix([1,-1,-1,0])
    check(J*gauge==sp.zeros(3,1),'Incorrect nuisance-null direction.')
    check(J.rank()==3,'Unexpected full-model rank.')
    check(J[:,[0,1,3]].rank()==3,'Restricted model should be identifiable.')
    return dict(contrast=str(expr),rank_full=int(J.rank()),rank_kappa_known=int(J[:,[0,1,3]].rank()),
                response_matrix=np.array(J).astype(float).tolist(),
                null_direction=[1,-1,-1,0],
                parameter_order=['beta','eta','kappa','lambda'],
                observation_order=['loop','initial_drift','first_order_spectrum'])


def bounded_beta_interval(loop: float, drift: float, loop_se: float, drift_se: float,
                          eta_bound: float, kappa_bound: float, alpha: float=.05,
                          drift_bias_bound: float=0.) -> dict[str,Any]:
    """Conditional >=1-alpha interval under asymptotic normal readout errors.

    Uses a Bonferroni rectangle for two channels and EXTERNALLY JUSTIFIED hard
    bounds |eta|<=eta_bound, |kappa|<=kappa_bound. Infinity means unconstrained.
    It does not infer or verify those bounds from the same experiment.
    Additional model/preparation/readout errors must be included, not ignored.
    """
    if not 0<alpha<1 or min(loop_se,drift_se,eta_bound,kappa_bound,drift_bias_bound)<0:
        raise ValueError('Invalid uncertainty bounds or significance.')
    z=float(norm.ppf(1-alpha/4))
    lo=max(loop-z*loop_se-eta_bound,
           drift-z*drift_se-kappa_bound-drift_bias_bound)
    hi=min(loop+z*loop_se+eta_bound,
           drift+z*drift_se+kappa_bound+drift_bias_bound)
    empty=lo>hi
    return dict(lower=None if not math.isfinite(lo) else lo,
                upper=None if not math.isfinite(hi) else hi,
                lower_unbounded=(lo==-math.inf),upper_unbounded=(hi==math.inf),
                empty=empty,zero_excluded=(not empty and (lo>0 or hi<0)),
                z=z,confidence_at_least=1-alpha,
                interpretation='Conditional on the chosen dynamics, asymptotic normal error model and independently valid nuisance/bias bounds. Empty means model/bounds inconsistency, not discovery.')


def run(c:Config, dest:Path) -> dict[str,Any]:
    check(c.dimension>=32 and c.repetitions>=2,'Insufficient dimension or repetitions.')
    check(c.fleet_trials>0 and c.fleet_trials%2==0,'Fleet trial count must be positive and even.')
    check(c.position_shots_per_setting>0 and 0<c.visibility<1,'Invalid shot count or visibility.')
    check(0<c.momenta[0]<c.momenta[1] and c.tmax>0 and c.time_points>=3,'Invalid controls.')
    dest.mkdir(parents=True,exist_ok=True)
    result:dict[str,Any]=dict(configuration=asdict(c),versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),algebra=algebra_checks())
    example=CASES[1]
    times=np.linspace(0,2,41)
    convergence=[]
    ref=trajectory(example,96,2.,times)
    for dim in (32,48,64,96):
        tr=trajectory(example,dim,2.,times)
        err=float(np.max(abs(tr['means']-ref['means'])))
        predicted=2.+4*example.beta/3*(2.**3+3*2.*.5)
        slope_error=abs(tr['initial_slope']-predicted)
        check(err<1e-8 and slope_error<1e-9,'Low-energy matrix convergence failed.')
        convergence.append(dict(dimension=dim,maximum_mean_error=err,initial_slope_error=slope_error,edge_probability=tr['edge_probability']))
    result['dimension_convergence']=convergence
    # Exact isospectral identity for the two quartic effective Hamiltonians.
    hp,x,p=hamiltonian(example,c.dimension)
    hx,_,_=hamiltonian(Case('spring',lam=example.beta),c.dimension)
    rotate=np.diag(np.exp(-.5j*np.pi*np.arange(c.dimension)))
    rotational_error=float(np.max(abs(rotate@hp@rotate.conj().T-hx)))
    ep=eigh(hp,eigvals_only=True);ex=eigh(hx,eigvals_only=True)
    spectral_error=float(np.max(abs(ep[:12]-ex[:12])))
    check(rotational_error<1e-10 and spectral_error<1e-10,'Quartic spectra should coincide.')
    result['isospectral_counterexample']=dict(rotation_matrix_error=rotational_error,
        maximum_low_12_eigenvalue_difference=spectral_error,
        kinetic_low_8_energies=ep[:8].tolist(),potential_low_8_energies=ex[:8].tolist(),
        first_order_spacing_slope=example.beta,
        exact_lowest_spacing_difference=float(ep[2]-2*ep[1]+ep[0]))
    # Strongest conventional mimic: identical free/pulse generators, not just moments.
    mimic=CASES[4]
    hm,_,_=hamiltonian(mimic,c.dimension)
    free_error=float(np.max(abs(hm-hp)))
    pulse_error=float(np.max(abs((example.loop-mimic.loop)*np.linalg.matrix_power(p,3)/3)))
    check(free_error==0 and pulse_error==0,'Complete mimic must have identical generators.')
    result['complete_mimic']=dict(free_hamiltonian_error=free_error,pulse_generator_error=pulse_error,
        statement='All protocols using these same preparations, measurement operators, free and control generators have identical distributions.')
    # Finite-window extrapolation: compare to the independent commutator derivative.
    windows=[]
    for case in [example,CASES[3],Case('deformation_plus_strong_spring',beta=.0025,lam=.01)]:
        for maximum in (.2,.4,.6,1.):
            a=drift_analysis(case,c,maximum)
            windows.append(dict(case=case.name,tmax=maximum,estimated_kinetic=a['estimate'],bias=a['bias'],predicted_se=a['predicted_se']))
    result['finite_window_audit']=windows
    # Known extra p^6 term tests the boundary of the first-order interpretation.
    h2,x,p=hamiltonian(example,c.dimension,True)
    ratios=[]
    for p0 in c.momenta:
        psi=packet(c.dimension,p0)
        ratios.append(float(np.vdot(psi,1j*(h2@x-x@h2)@psi).real)/p0)
    second_drift=3*(ratios[1]-ratios[0])/(4*(c.momenta[1]**2-c.momenta[0]**2))
    second_expected=example.beta+17/20*example.beta**2*(sum(q*q for q in c.momenta)+5)
    check(abs(second_drift-second_expected)<1e-10,'Second-order derivative calculation failed.')
    result['higher_order_audit']=dict(first_order=example.beta,with_p6=second_drift,
                                    fractional_difference=second_drift/example.beta-1,
                                    caveat='Retains the next local tangent expansion term, not a global minimal-length operator theory.')
    # A nonlinear kick calibration can mimic the initial-slope contrast alone.
    zeta=4*example.beta/3
    kick_mimic=3*zeta/4
    result['preparation_calibration_counterexample']=dict(actual_mean_momentum='q+zeta*q^3',
        zeta=zeta,inferred_kinetic_if_nominal_q_is_used=kick_mimic,
        fractional_momentum_error_at_q1=zeta*c.momenta[0]**2,
        fractional_momentum_error_at_q2=zeta*c.momenta[1]**2,
        scope='Mimics the slope contrast; not a claim of identical full trajectories.')
    rng=np.random.default_rng(c.seed)
    samples=[]
    for case in CASES:
        print('Sampling',case.name,flush=True)
        sample=sample_case(case,c,rng)
        samples.append(sample)
        print('  drift',sample['estimated_kinetic']['mean'],'+/-',sample['estimated_kinetic']['standard_deviation'],flush=True)
    result['sampling']=samples
    result['resources_per_simulated_experiment']=dict(
        position_settings=2*len(c.momenta)*c.time_points,
        destructive_position_readings=2*len(c.momenta)*c.time_points*c.position_shots_per_setting,
        additional_ghz_preparations=c.fleet_trials,additional_local_qubit_readings=3*c.fleet_trials,
        unit_time='tau=omega*t',
        note='No equal-resource sensitivity advantage is claimed. Calibration, state preparation, readout engineering and dissipative costs are not modeled.')
    # Conditional calibration floors, in the same dimensionless coefficient units.
    sd=samples[1]['estimated_kinetic']['analytic_standard_error']
    kappas=[0.,.0001,.001,.01]
    result['independent_calibration_floors']=[dict(independent_kappa_se=s,
        beta_standard_error=math.hypot(sd,s),
        infinite_shot_floor=s) for s in kappas]
    loop_se=samples[1]['estimated_loop']['analytic_standard_error']
    result['calibration_certificate_examples']=[]
    for label,eb,kb in [('no_independent_bounds',math.inf,math.inf),
                        ('assumed_loose_bounds',.003,.003),
                        ('assumed_tight_bounds',.0001,.0001)]:
        certificate=bounded_beta_interval(.0025,.0025,loop_se,sd,eb,kb,
                                         drift_bias_bound=.000001)
        result['calibration_certificate_examples'].append(dict(
            label=label,eta_bound=None if math.isinf(eb) else eb,
            kappa_bound=None if math.isinf(kb) else kb,
            bounds_are_hypothetical_not_experimental=True,**certificate))
    result['status']=dict(algebra_passed=True,matrix_checks_passed=True,
        actuator_only_distinguishable_under_completion=True,
        nonlinear_spring_distinguishable_by_initial_drift=True,
        complete_conventional_mimic_distinguishable=False,
        physically_validated_subplanck_sensitivity=False)
    (dest/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (dest/'frozen_parameters.json').write_text(json.dumps(asdict(c),indent=2)+'\n')
    with (dest/'recovery_summary.csv').open('w',newline='') as f:
        names=['case','true_beta','true_eta','true_kappa','true_lambda','estimated_loop','loop_sd','estimated_kinetic','kinetic_sd','conditional_eta','conditional_eta_sd']
        writer=csv.DictWriter(f,fieldnames=names);writer.writeheader()
        for s in samples:
            writer.writerow(dict(case=s['case']['name'],true_beta=s['case']['beta'],true_eta=s['case']['eta'],
                true_kappa=s['case']['kappa'],true_lambda=s['case']['lam'],
                estimated_loop=s['estimated_loop']['mean'],loop_sd=s['estimated_loop']['standard_deviation'],
                estimated_kinetic=s['estimated_kinetic']['mean'],kinetic_sd=s['estimated_kinetic']['standard_deviation'],
                conditional_eta=s['conditional_eta']['mean'],conditional_eta_sd=s['conditional_eta']['standard_deviation']))
    with (dest/'trajectories.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['case','momentum','tau','expected_position','position_variance'])
        for case in CASES:
            for p0 in c.momenta:
                tr=trajectory(case,c.dimension,p0,times)
                for t,m,v in zip(times,tr['means'],tr['variances']):
                    writer.writerow([case.name,p0,t,m,v])
    print(json.dumps(result['status'],indent=2),flush=True)
    return result


def main()->None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('.'))
    parser.add_argument('--repetitions',type=int,default=5000)
    parser.add_argument('--position-shots',type=int,default=1_000_000)
    parser.add_argument('--seed',type=int,default=20260906)
    args=parser.parse_args()
    run(Config(repetitions=args.repetitions,position_shots_per_setting=args.position_shots,seed=args.seed),args.output)

if __name__=='__main__':
    main()
