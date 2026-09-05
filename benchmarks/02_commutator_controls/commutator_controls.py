#!/usr/bin/env python3
"""Reproducible ideal-control sensor commutator benchmark, v0.2.

This is a conditional quantum-model calculation, NOT a hardware experiment,
not a UV completion, and not an observation of a minimal length.

Dependencies: numpy, scipy, sympy. All pulses are dimensionless integrated
Hamiltonian strengths. The first-order GUP is represented by the ordinary
canonical polynomial P=p+(beta+eta)*p**3/3. eta is an indistinguishable
conventional actuator nonlinearity. The polynomial model itself is exact;
its identification with the specified GUP is only first order.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
import math
from dataclasses import dataclass, asdict, replace
from pathlib import Path
from typing import Any
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.special import roots_hermitenorm
import sympy as sp

SIGNS = np.array([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]], dtype=int)
ZSTATES = np.array(list(itertools.product([1,-1], repeat=3)), dtype=int)

@dataclass(frozen=True)
class Config:
    a: float = 0.5
    b: float = 0.5
    q0: float = 1.0
    g: tuple[float,float,float] = (1.0,1.0,1.0)
    beta: float = 0.0025
    eta: float = 0.0
    visibility: float = 0.8
    trials_per_experiment: int = 100000
    repetitions: int = 20000
    seed: int = 20260905


def leading_phase(c: Config) -> float:
    return float(128*(c.beta+c.eta)*c.a**3*c.b*c.q0*np.prod(c.g))


def symmetric_loop_phase(p: np.ndarray, q: float, a: float, b: float,
                         beta: float, eta: float=0.0, fifth: bool=False) -> np.ndarray:
    """Exact phase of the signed loop pair for the chosen momentum polynomial."""
    be=beta+eta
    out=-2*a*b*q*q*(1+be*p*p)-(2/3)*be*a**3*b*q**4
    if fifth:
        # P=p+beta*p^3/3+2*beta^2*p^5/15: second-order tan-series test.
        if eta != 0:
            raise ValueError('Fifth-order check uses beta alone.')
        out += -(4/3)*beta**2*a*b*q*q*p**4
        out += -(8/3)*beta**2*a**3*b*q**4*p*p
        out += -(4/15)*beta**2*a**5*b*q**6
    return out


def branch_phase(c: Config, z: np.ndarray, p: np.ndarray,
                 imbalance: float=0., fifth: bool=False,
                 nuisance: np.ndarray | None=None) -> np.ndarray:
    phase=np.zeros_like(p, dtype=float)
    for k,s in enumerate(SIGNS):
        q=float(c.q0+np.dot(np.asarray(c.g),s*z))
        aa=c.a*(1+imbalance if k==0 else 1.)
        phase += symmetric_loop_phase(p,q,aa,c.b,c.beta,c.eta,fifth)
        if nuisance is not None:
            sz=s*z
            vals=np.array([1,sz[0],sz[1],sz[2],sz[0]*sz[1],sz[0]*sz[2],sz[1]*sz[2]])
            phase -= float(np.dot(nuisance,vals))
    return phase


def gaussian_coherence(c: Config, variance: float=0.5,
                       imbalance: float=0., fifth: bool=False,
                       nuisance: np.ndarray | None=None) -> complex:
    """Trace oscillator for a centered Gaussian momentum marginal."""
    if variance<=0:
        raise ValueError('Momentum variance must be positive.')
    nodes,w=roots_hermitenorm(96)
    p=nodes*np.sqrt(variance)
    pp=branch_phase(c,np.ones(3,int),p,imbalance,fifth,nuisance)
    pm=branch_phase(c,-np.ones(3,int),p,imbalance,fifth,nuisance)
    return complex(np.dot(w,np.exp(1j*(pm-pp)))/np.sqrt(2*np.pi))


def matrix_coherence(c: Config, dim: int=64, occupation: float=0.0) -> complex:
    """Independent multiplication of all 32 oscillator pulses, not phase formula.

    A truncated Fock basis cannot satisfy an exact canonical commutator globally.
    Dimension convergence is therefore checked on the specified low-energy state.
    """
    an=np.diag(np.sqrt(np.arange(1,dim)),1)
    x=(an+an.T)/np.sqrt(2)
    pc=(an-an.T)/(1j*np.sqrt(2))
    ex,vx=eigh(x); ep,vp=eigh(pc)
    pp=ep+(c.beta+c.eta)*ep**3/3
    def exp_i(vals: np.ndarray, v: np.ndarray, t: float) -> np.ndarray:
        return (v*np.exp(1j*t*vals))@v.conj().T
    def loop(a: float,b: float,q: float) -> np.ndarray:
        # Rightmost pulse acts first.
        return exp_i(ex,vx,a*q)@exp_i(pp,vp,b*q)@exp_i(ex,vx,-a*q)@exp_i(pp,vp,-b*q)
    cache: dict[float,np.ndarray]={}
    def sym(q: float) -> np.ndarray:
        if q not in cache:
            cache[q]=loop(-c.a,-c.b,q)@loop(c.a,c.b,q)
        return cache[q]
    evol=[]
    for z in (np.ones(3,int),-np.ones(3,int)):
        u=np.eye(dim,dtype=complex)
        for s in SIGNS:
            q=float(c.q0+np.dot(c.g,s*z))
            u=sym(q)@u
        evol.append(u)
    if occupation==0:
        weights=np.zeros(dim); weights[0]=1
    else:
        r=occupation/(occupation+1.)
        weights=(1-r)*r**np.arange(dim); weights/=sum(weights)
    return complex(np.sum(weights*np.sum(evol[0].conj()*evol[1],axis=0)))


def symbolic_check() -> dict[str,str]:
    q0,g1,g2,g3,a,b,p,be=sp.symbols('q0 g1 g2 g3 a b p beta',real=True)
    branch=[]
    for z in ZSTATES:
        val=0
        for s in SIGNS:
            q=q0+sum(gi*int(si)*int(zi) for gi,si,zi in zip([g1,g2,g3],s,z))
            val+=-2*a*b*q**2*(1+be*p*p)-sp.Rational(2,3)*be*a**3*b*q**4
        branch.append(sp.expand(val))
    result={}
    for mask in itertools.product([0,1],repeat=3):
        coeff=sp.factor(sum(e*int(np.prod(z**np.array(mask))) for z,e in zip(ZSTATES,branch))/8)
        name=''.join(map(str,mask)); result[name]=str(coeff)
        if name not in ('000','111'):
            assert coeff==0, (name,coeff)
        if name=='111':
            assert sp.expand(coeff+64*a**3*b*be*g1*g2*g3*q0)==0
    return result


def run(c: Config, dest: Path) -> dict[str,Any]:
    if not 0<c.visibility<=1 or c.trials_per_experiment<2 or c.repetitions<2:
        raise ValueError('Invalid visibility or sample counts.')
    if c.trials_per_experiment%2:
        raise ValueError('Total trials must be even for the two-setting experiment.')
    dest.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(c.seed)
    sym=symbolic_check()
    matrix=[]
    for dim in (32,48,64,96):
        cc=matrix_coherence(c,dim)
        err=abs(cc-np.exp(1j*leading_phase(c)))
        assert err<1e-9
        matrix.append(dict(dimension=dim,phase=float(np.angle(cc)),
                           coherence_modulus=float(abs(cc)),complex_error=float(err)))
    general=replace(c,a=.45,b=.3,q0=.7,g=(.4,.8,1.1),beta=.004)
    gc=matrix_coherence(general,64)
    assert abs(gc-np.exp(1j*leading_phase(general)))<1e-9
    thermal=matrix_coherence(c,96,occupation=1.)
    assert abs(thermal-np.exp(1j*leading_phase(c)))<1e-8
    pure=replace(c,beta=0.,eta=0.)
    scenarios=[('canonical_ideal',pure,0.),
               ('specified_deformation',c,0.),
               ('ordinary_cubic_actuator_mimic',replace(c,beta=0.,eta=c.beta),0.),
               ('bias_reversed',replace(c,q0=-c.q0),0.),
               ('bias_off',replace(c,q0=0.),0.),
               ('one_sensor_decoupled',replace(c,g=(1.,1.,0.)),0.),
               ('canonical_one_block_0.1pct_imbalance',pure,.001)]
    rows=[]
    for label,case,eps in scenarios:
        coh=gaussian_coherence(case,imbalance=eps)
        rows.append(dict(case=label,phase_radians=float(np.angle(coh)),
                         coherence_modulus=float(abs(coh)),
                         expected_parity=float(c.visibility*coh.imag)))
    assert abs(rows[1]['phase_radians']-rows[2]['phase_radians'])<1e-13
    assert abs(rows[-1]['phase_radians']-.006)<1e-12
    nuisance_errors=[]
    for _ in range(100):
        nuisance=rng.uniform(-2.,2.,size=7)
        coh=gaussian_coherence(c,nuisance=nuisance)
        nuisance_errors.append(abs(coh-np.exp(1j*leading_phase(c))))
    assert max(nuisance_errors)<1e-11
    # Second-order model-error test, not claimed to be an all-orders GUP theory.
    higher=[]
    for variance in (.5,2.5,10.5):
        cc=gaussian_coherence(c,variance=variance,fifth=True)
        expected_phase=leading_phase(c)+1024*c.beta**2*c.a**5*c.b+0.5*np.arctan(1024*c.beta**2*c.a**3*c.b*variance)
        expected_modulus=(1+(1024*c.beta**2*c.a**3*c.b*variance)**2)**(-.25)
        assert abs(np.angle(cc)-expected_phase)<1e-11
        assert abs(abs(cc)-expected_modulus)<1e-11
        higher.append(dict(momentum_variance=variance,phase=float(np.angle(cc)),
                           coherence_modulus=float(abs(cc)),relative_phase_change=float(np.angle(cc)/leading_phase(c)-1)))
    # Sampling the derived probability law; bias reversal is a physical q0 change.
    # Fixed common readout offset is removed, but the actuator mimic reverses too.
    offset=.03
    n=c.trials_per_experiment//2
    monte=[]
    for label,case in [('deformation',c),('actuator_mimic',replace(c,beta=0.,eta=c.beta)),('null',pure)]:
        phi=leading_phase(case)
        kp=rng.binomial(n,(1+c.visibility*np.sin(offset+phi))/2,size=c.repetitions)
        km=rng.binomial(n,(1+c.visibility*np.sin(offset-phi))/2,size=c.repetitions)
        hp=np.arcsin(np.clip((2*kp/n-1)/c.visibility,-1,1))
        hm=np.arcsin(np.clip((2*km/n-1)/c.visibility,-1,1))
        phi_est=(hp-hm)/2
        slope=128*c.a**3*c.b*c.q0*np.prod(c.g)
        est=phi_est/slope
        monte.append(dict(case=label,true_effective_beta=float(case.beta+case.eta),
                          mean_effective_beta=float(est.mean()),std_effective_beta=float(est.std(ddof=1)),
                          mean_recovered_phase=float(phi_est.mean()),
                          total_triples_per_experiment=c.trials_per_experiment,
                          repetitions=c.repetitions))
    # Full recorded triples: exact parity law, uniform proper marginals.
    raw_n=c.trials_per_experiment
    phi=leading_phase(c)
    w=np.where(rng.random(raw_n)<(1+c.visibility*np.sin(phi))/2,1,-1)
    x1=2*rng.integers(0,2,raw_n)-1; x2=2*rng.integers(0,2,raw_n)-1
    x3=w*x1*x2
    records=np.column_stack([x1,x2,x3])
    with (dest/'raw_triple_counts.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['x1','x2','x3','count'])
        for z in ZSTATES:writer.writerow([*z,int(np.sum(np.all(records==z,axis=1)))])
    # Multiple controls identify effective beta vs simple imbalance but not beta vs eta.
    controls=[]
    for u in (.7,1.,1.3):
        for s in (-1,1):
            controls.append([s*8*u**4,s*8*u**4,s*6*u**2,1.])
    design=np.array(controls)
    singular=np.linalg.svd(design,compute_uv=False)
    rank=int(np.linalg.matrix_rank(design))
    assert rank==3
    assert np.linalg.matrix_rank(design[:,[0,2,3]])==3
    # Unit conversion. Assumed oscillator values, NOT measured apparatus parameters.
    h=6.62607015e-34; hbar=h/(2*np.pi); light=299792458.
    mass_planck=2.176434e-8; length_planck=1.616255e-35
    mass=1e-11; freq=1e5; omega=2*np.pi*freq
    pstar=float(np.sqrt(hbar*mass*omega)); xstar=float(np.sqrt(hbar/(mass*omega)))
    ratio=float((pstar/(mass_planck*light))**2)
    slope=128*c.a**3*c.b*c.q0*np.prod(c.g)
    scale=[]
    for beta0 in (1.,.01):
        beta=float(beta0*ratio); phi=float(slope*beta)
        n5=float(25/(c.visibility**2*phi**2))
        scale.append(dict(beta0=beta0,model_minimal_length_m=float(np.sqrt(beta0)*length_planck),
                          dimensionless_beta=beta,phase_radians=phi,
                          expected_5_null_standard_error_trials=n5))
    target=.02
    # a=b=A, q0=g_i=1. Raw equal-area cancellation tolerance in one block.
    A=float((target/(128*ratio))**.25)
    tolerance=float(target/(24*A*A))
    output=dict(version='0.2',status='Conditional derivation and simulation; no hardware or sub-Planck observation',
                config=asdict(c),versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
                symbolic_walsh_coefficients=sym,matrix_convergence=matrix,
                unequal_coupling_check=dict(config=asdict(general),expected_phase=leading_phase(general),
                                            matrix_phase=float(np.angle(gc)),error=float(abs(gc-np.exp(1j*leading_phase(general))))),
                thermal_matrix_check=dict(mean_occupation=1.,dimension=96,phase=float(np.angle(thermal)),
                                          error=float(abs(thermal-np.exp(1j*leading_phase(c))))),
                cases=rows,max_static_diagonal_nuisance_error=float(max(nuisance_errors)),
                second_order_momentum_model=higher,monte_carlo=monte,
                raw_record_summary=dict(triples=raw_n,single_means=records.mean(axis=0).tolist(),
                    pair_means=[float((x1*x2).mean()),float((x1*x3).mean()),float((x2*x3).mean())],
                    parity_mean=float(w.mean())),
                identifiability=dict(parameter_order=['beta','eta','block_imbalance','common_readout_offset'],
                    rank=rank,parameter_count=4,singular_values=singular.tolist(),
                    reduced_effective_parameter_rank=3,
                    interpretation='All measured controls depend on beta+eta, never beta and eta separately.'),
                units=dict(assumed_mass_kg=mass,assumed_frequency_Hz=freq,pstar_kg_m_per_s=pstar,
                    xstar_m=xstar,planck_mass_kg=mass_planck,planck_length_m=length_planck,
                    beta_per_unit_beta0=ratio,toy_beta0=float(c.beta/ratio)),
                length_scale_examples=scale,
                amplification_illustration=dict(target_phase=target,a_equals_b=A,
                    largest_conditional_kick=4*A,corresponding_position_translation_m=4*A*xstar,
                    one_block_fractional_a_error_for_signal_sized_leakage=tolerance,
                    status='Ideal algebra only; no feasible implementation or precision demonstrated.'),
                resources=dict(oscillator_pulses_per_trial=32,symmetric_loop_pairs=4,
                    individual_qubit_X_gates_for_closed_toggling_cycle=8,
                    preparations_per_two_setting_experiment=c.trials_per_experiment,
                    local_readings_per_two_setting_experiment=3*c.trials_per_experiment,
                    oscillator_pulses_per_two_setting_experiment=32*c.trials_per_experiment),
                checks_passed=True)
    (dest/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    with (dest/'scenario_results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return output


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--seed',type=int,default=20260905)
    parser.add_argument('--repetitions',type=int,default=20000)
    args=parser.parse_args()
    result=run(Config(seed=args.seed,repetitions=args.repetitions),args.output)
    print(json.dumps({k:result[k] for k in ['checks_passed','cases','monte_carlo','length_scale_examples','amplification_illustration']},indent=2))

if __name__=='__main__':
    main()
