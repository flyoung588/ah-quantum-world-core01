"""Actual one-electron Coulomb eigenproblem and midpoint electronic dipole.
No SCF/FCI correlation solver. Sigma identity is checked with Lz action, not energy order alone.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,solve

def invert_AO(mol):
    n=mol.nao_nr();assert n%2==0;loc=mol.ao_loc_nr();parity=np.empty(n)
    for sh in range(mol.nbas):parity[loc[sh]:loc[sh+1]]=(-1)**mol.bas_angular(sh)
    h=n//2;assert np.array_equal(parity[:h],parity[h:]);P=np.zeros((n,n));P[:h,h:]=np.diag(parity[:h]);P[h:,:h]=np.diag(parity[:h]);return P

def point(R,basis,previous=None):
    from pyscf import gto
    mol=gto.M(atom=[['H',(0.,0.,-R/2)],['H',(0.,0.,R/2)]],basis=basis,unit='Bohr',charge=1,spin=1,cart=False,verbose=0,max_memory=400)
    assert mol.nelectron==1 and mol.natm==2 and abs(mol.energy_nuc()-1/R)<1e-12
    mol.set_common_origin((0.,0.,0.))
    S=mol.intor('int1e_ovlp');h=mol.intor('int1e_kin')+mol.intor('int1e_nuc');dip=-mol.intor('int1e_r',comp=3)
    Lz=mol.intor('int1e_cg_irxp',comp=3)[2];P=invert_AO(mol)
    eigS=np.linalg.eigvalsh(S);assert eigS[0]>1e-10,'AO near linear dependence: no silent truncation'
    eps,C=eigh(h,S);parities=np.sum(C*(S@P@C),axis=0)
    action=solve(S,Lz@C,assume_a='pos');sigma=np.sqrt(np.maximum(np.sum(action*(S@action),axis=0),0))
    indices=[]
    for target in [1,-1]:
        candidates=np.flatnonzero((abs(parities-target)<1e-7)&(sigma<1e-7))
        if not len(candidates):raise RuntimeError('no symmetry-identified sigma root')
        indices.append(int(candidates[0]))
    coeff=C[:,indices].copy();signs=[1.,1.];overlaps=[None,None]
    if previous is not None:
        oldmol,oldcoeff=previous;cross=gto.intor_cross('int1e_ovlp',oldmol,mol)
        for k in range(2):
            overlap=float(oldcoeff[:,k]@cross@coeff[:,k]);signs[k]=1. if overlap>=0 else -1.;coeff[:,k]*=signs[k];overlaps[k]=abs(overlap)
            if abs(overlap)<.7:raise RuntimeError('electronic phase/root continuation not trustworthy')
    else:
        c0=(mol.eval_gto('GTOval_sph',np.zeros((1,3)))@coeff[:,0]).item();signs[0]=1. if c0>=0 else -1.;coeff[:,0]*=signs[0]
        d=float(coeff[:,0]@dip[2]@coeff[:,1]);signs[1]=1. if d>=0 else -1.;coeff[:,1]*=signs[1]
    transitions=np.array([coeff[:,0]@x@coeff[:,1] for x in dip]);res=[]
    for k,j in enumerate(indices):
        error=h@coeff[:,k]-eps[j]*S@coeff[:,k];res.append(float(np.sqrt(max(error@solve(S,error,assume_a='pos'),0))))
    assert max(res)<1e-9 and max(abs(transitions[:2]))<1e-8
    row={'R_bohr':float(R),'basis':basis,'AO_count':mol.nao_nr(),'AO_overlap_min_eigenvalue':float(eigS[0]),'electronic_energy_Eh':eps[indices].tolist(),'nuclear_repulsion_Eh':1/float(R),'Vg_Eh':float(eps[indices[0]]+1/R),'Vu_Eh':float(eps[indices[1]]+1/R),'mu_gu_au':float(transitions[2]),'dipole_xyz_au':transitions.tolist(),'dipole_convention':'electronic -r, nuclear midpoint origin; nuclei transition term zero','root_indices':indices,'inversion_parity':parities[indices].tolist(),'Lz_action_norm':sigma[indices].tolist(),'residual_S_inverse_norm_Eh':res,'phase_sign_applied':signs,'adjacent_phase_overlap':overlaps,'AO_inversion_metric_error':float(np.max(abs(P.T@S@P-S))),'normalization_CtSC':(coeff.T@S@coeff).tolist()}
    return row,(mol,coeff),{'S':S,'h':h,'dipole':dip,'Lz':Lz,'coefficients':coeff,'AO_labels':mol.ao_labels()}

def build_curve(basis,directory):
    p=Path(directory);p.mkdir(exist_ok=True,parents=True);rows=[];previous=None;coefficients=[]
    for R in np.round(np.arange(.7,12.0001,.05),8):
        row,previous,raw=point(float(R),basis,previous);rows.append(row);coefficients.append(raw['coefficients'])
        if any(abs(R-x)<1e-8 for x in [1.,2.,3.,4.,6.,8.,10.]):np.savez_compressed(p/('AO_R'+str(R)+'.npz'),S=raw['S'],h=raw['h'],dipole=raw['dipole'],Lz=raw['Lz'],coefficients=raw['coefficients'])
        tmp=p/'curve.partial.json';tmp.write_text(json.dumps({'basis':basis,'rows':rows},indent=2))
    data={'basis':basis,'source':'PySCF one-electron AO Coulomb integrals, generalized Hermitian eigh, inversion and axial-angular identity checks','states':['1s sigma_g','2p sigma_u'],'origin_bohr':[0,0,0],'rows':rows}
    (p/'curve.json').write_text(json.dumps(data,indent=2));np.savez_compressed(p/'electronic_coefficients.npz',coefficients=np.array(coefficients))
    return data
