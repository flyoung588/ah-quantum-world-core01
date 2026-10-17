"""Tiny synthetic arithmetic/branch/checkpoint tests; NEVER actual H2+ data."""
import json,sys,tempfile,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'core02'))
import numpy as np
from runtime import QuantumWorld,PhysicalSystem,QuantumState
from coupled import MatrixHamiltonian,local_exponential,field
start=time.monotonic();n=7;R=np.linspace(1.,2.,n)
class Toy(PhysicalSystem):
 max_dt=.1
 def __init__(self):
  self.axes=[R];self.dv=R[1]-R[0];self.hamiltonian=MatrixHamiltonian(R,np.full(n,.1),np.full(n,.3),np.ones(n),1000.)
 def metadata(self):return {'synthetic_software_fixture':True,'n':n}
 def allowed(self):return {'toy':True}
 def initial_wavefunction(self,v):return np.stack([np.sin(np.pi*np.arange(1,n+1)/(n+1)),np.zeros(n)])
 def apply_to_state(self,s,c):s.parameters['pulse']=dict(start_au=s.t,amplitude_au=c['amplitude'],duration_au=50.,omega_au=.2,cep_rad=0.)
 def observables(self,s):return {'Pu':float(np.sum(abs(s.psi[1])**2)*self.dv)}
psi=np.array([[1.+0j],[0j]]);coupling=.03;detuning=.2;timepoint=3.;value=local_exponential(psi,np.array([.1]),np.array([.3]),np.array([coupling]),timepoint)
omega=np.sqrt((detuning/2)**2+coupling**2);expected=coupling**2/omega**2*np.sin(omega*timepoint)**2
error=abs(abs(value[1,0])**2-expected);assert error<1e-14
w=QuantumWorld(Toy());w.prepare({});branch=w.fork();branch.apply({'amplitude':.04})
assert w.state.parameters=={} and branch.state.parameters!=w.state.parameters
w.evolve(1.,.1);branch.evolve(1.,.1);assert w.observe()['Pu']==0 and branch.observe()['Pu']>0
assert abs(branch.observe()['norm']-1)<1e-12
with tempfile.TemporaryDirectory(prefix='core03-synthetic-') as tmp:
 branch.checkpoint(tmp);restored=QuantumWorld(branch.system);restored.restore(tmp)
 assert np.array_equal(restored.state.psi,branch.state.psi) and restored.state.parameters==branch.state.parameters and restored.state.work_input==branch.state.work_input
pulse=branch.state.parameters['pulse'];t=13.;eps=1e-5;derivative_error=abs(field(t,pulse)[1]-(field(t+eps,pulse)[0]-field(t-eps,pulse)[0])/(2*eps));assert derivative_error<1e-10
assert 'pyscf' not in sys.modules
result={'scope':'SYNTHETIC_SOFTWARE_ONLY, seven radial grid points; no molecular physical computation','Rabi_formula_error':float(error),'norm_error':abs(branch.observe()['norm']-1),'state_owned_control_isolation':'PASS','checkpoint_restores_state_control_and_work':'PASS','field_derivative_error':derivative_error,'wall_seconds':time.monotonic()-start,'physical_accuracy':'NOT_TESTED','PySCF_imported':False}
print(json.dumps(result,indent=2))
