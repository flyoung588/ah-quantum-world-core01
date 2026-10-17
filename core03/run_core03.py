"""Single bounded cloud experiment, after explicit external authorization only."""
import hashlib,json,os,platform,resource,sys,time,traceback,subprocess
from pathlib import Path
P=Path(__file__).parent;sys.path[:0]=[str(P),str(P.parent/'core02'),str(P.parent)]
import numpy as np,scipy
from runtime import QuantumWorld,clean
from electronic import build_curve
from h2plus_adapter import H2PlusSystem
from coupled import local_exponential
OUT=Path('/out');OUT.mkdir(exist_ok=True)
def save(name,value):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(clean(value),indent=2,allow_nan=False));tmp.replace(p)
def cgroup():return {k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in ['memory.max','memory.current','memory.peak','memory.events','memory.swap.max','cpu.max','cpu.stat']}
def propagate(system,initial,pulse,duration,dt,name=None,prepared_world=None):
 w=QuantumWorld(system) if prepared_world is None else prepared_world
 if w.state is None:w.prepare(initial)
 if 'pulse' not in w.state.parameters:w.apply(pulse)
 base=w.observe();rows=[base];density=[]
 for i in range(round(duration/(10*dt))):
  w.step(dt,count=10);o=w.observe();rows.append(o)
  if i%12==0:density.append({'t_au':w.state.t,'g':(abs(w.state.psi[0])**2).tolist(),'u':(abs(w.state.psi[1])**2).tolist()})
 result={'initial':initial,'resolved_pulse':w.state.parameters['pulse'],'rows':rows,'density_frames':density,'energy_change_Eh':rows[-1]['energy_Eh']-base['energy_Eh'],'external_work_Eh':w.state.work_input,'work_energy_error_Eh':rows[-1]['energy_Eh']-base['energy_Eh']-w.state.work_input,'total_probability_max_error':max(abs(o['Pg']+o['Pu']-1) for o in rows),'max_edge_probability':max(o['edge_probability'] for o in rows),'classification':'PROGRAMMATIC_EXPERIMENT'}
 if name:
  save(name+'.json',result);w.error_evidence={'probability_max_error':result['total_probability_max_error'],'work_energy_error_Eh':result['work_energy_error_Eh'],'reference':'sensitivity.json; not a physical-accuracy certificate'};w.export_evidence(OUT/(name+'_state'))
 return w,result
def main():
 start=time.monotonic();cg=cgroup();assert cg['memory.max']=='1073741824' and cg['memory.swap.max']=='0' and cg['cpu.max']=='100000 100000'
 assert os.listdir('/sys/class/net')==['lo'] and len(os.sched_getaffinity(0))==1
 save('admission.json',{'cgroup':cg,'network':os.listdir('/sys/class/net'),'affinity':list(os.sched_getaffinity(0))})
 import pyscf
 pyscf.lib.num_threads(1)
 save('environment.json',{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pyscf':pyscf.__version__,'sha':os.environ.get('GITHUB_SHA'),'pyscf_threads':pyscf.lib.num_threads()})
 cfg=json.loads((P/'requests.json').read_text());save('actual_input.json',cfg)
 for basis in [cfg['primary_basis'],cfg['secondary_basis']]:build_curve(basis,OUT/basis)
 primary=OUT/cfg['primary_basis']/'curve.json';secondary=OUT/cfg['secondary_basis']/'curve.json'
 system=H2PlusSystem(primary,cfg['nuclear_grid']);initial=cfg['initial'];pulse=cfg['pulse'];duration=cfg['duration_au'];dt=cfg['dt_au']
 center=initial['center_R_bohr'];pulse={**pulse,'omega_au':float(system.curves['Vu_Eh'](center)-system.curves['Vg_Eh'](center))*pulse['omega_ratio_to_local_gap']}
 save('resolved_control.json',{'pulse':pulse,'frequency_held_fixed_for_grid_interpolation_and_basis_checks':True})
 zero={**pulse,'amplitude_au':0.};A,a=propagate(system,initial,zero,duration,dt,'A_zero_field')
 # Fork from the SAME actual prepared wavefunction, protocol remains state-owned.
 root=QuantumWorld(system);root.prepare(initial);B=root.fork();C=root.fork();B.apply(pulse);branchpulse={**pulse,'amplitude_au':cfg['branch_amplitude_au']};C.apply(branchpulse)
 assert root.state.parameters=={} and B.state.parameters['pulse']['amplitude_au']!=C.state.parameters['pulse']['amplitude_au']
 # Store complete traces using the same generic interface, not separate model loops.
 B,b=propagate(system,initial,pulse,duration,dt,'B_driven',B);C,c=propagate(system,initial,branchpulse,duration,dt,'C_amplitude_branch',C)
 save('branch_comparison.json',{'comparison':B.compare(C),'same_prepared_source':True,'changed_parameter_only':'amplitude_au','zero_field_final':A.observe(),'no_shared_mutable_pulse':root.state.parameters=={}})
 assert max(o['Pu'] for o in a['rows'])<1e-14 and max(o['Pu'] for o in b['rows'])>1e-6
 assert max(v['total_probability_max_error'] for v in [a,b,c])<1e-9
 assert abs(a['energy_change_Eh'])<1e-6
 # Fixed-R constant physical electric coupling: exact detuned two-level Rabi formula, no RWA claim.
 R=2.;vg=float(system.curves['Vg_Eh'](R));vu=float(system.curves['Vu_Eh'](R));mu=float(system.curves['mu_gu_au'](R));E=.004;t=123.
 state=local_exponential(np.array([[1.+0j],[0j]]),np.array([vg]),np.array([vu]),np.array([-mu*E]),t)
 from scipy.linalg import expm
 matrix=np.array([[vg,-mu*E],[-mu*E,vu]]);alternate=expm(-1j*matrix*t)@np.array([1.,0.])
 omega=np.sqrt(((vg-vu)/2)**2+(mu*E)**2);expected=(mu*E)**2/omega**2*np.sin(omega*t)**2
 save('fixed_R_Rabi.json',{'R_bohr':R,'physical_mu_gu_au':mu,'Vg_Eh':vg,'Vu_Eh':vu,'constant_field_au':E,'time_au':t,'computed_Pu':float(abs(state[1,0])**2),'analytic_Pu':float(expected),'absolute_error':float(abs(abs(state[1,0])**2-expected)),'alternate_expm_max_error':float(np.max(abs(state[:,0]-alternate))),'Hermiticity_error':float(np.max(abs(matrix-matrix.T.conj()))),'meaning':'exact detuned constant-H two-level oscillation; not resonant-CW RWA, not independent dipole validation'})
 assert abs(abs(state[1,0])**2-expected)<1e-12
 dtworld,d=propagate(system,initial,pulse,duration,dt/2)
 gridworld,g=propagate(H2PlusSystem(primary,cfg['refined_nuclear_grid']),initial,pulse,duration,dt)
 interpworld,i=propagate(H2PlusSystem(primary,cfg['nuclear_grid'],stride=2),initial,pulse,duration,dt)
 basisworld,q=propagate(H2PlusSystem(secondary,cfg['nuclear_grid']),initial,pulse,duration,dt)
 def delta(result):return {'final_Pu_difference':result['rows'][-1]['Pu']-b['rows'][-1]['Pu'],'mean_R_difference':result['rows'][-1]['mean_R']-b['rows'][-1]['mean_R'],'work_energy_error_Eh':result['work_energy_error_Eh']}
 primary_rows=json.loads(primary.read_text())['rows'];secondary_rows=json.loads(secondary.read_text())['rows'];electronic={}
 for key in ['Vg_Eh','Vu_Eh','mu_gu_au']:
  differences=np.array([y[key]-x[key] for x,y in zip(primary_rows,secondary_rows)])
  maximum=int(np.argmax(abs(differences)))
  held_out=primary_rows[1::2];coarse_interp=interpworld.system.curves[key]
  interpolation_error=max(abs(float(coarse_interp(x['R_bohr']))-x[key]) for x in held_out)
  electronic[key]={'larger_basis_max_absolute_difference':float(max(abs(differences))),'at_R_bohr':primary_rows[maximum]['R_bohr'],'coarse_spline_vs_held_out_actual_AO_points_max_error':interpolation_error}
 save('sensitivity.json',{'time_step_half':delta(d),'radial_grid_doubled':delta(g),'electronic_interpolation_coarsened':delta(i),'larger_basis':delta(q),'actual_electronic_potential_dipole_checks':electronic,'baseline_work_energy_error_Eh':b['work_energy_error_Eh'],'basis_interpolation_differences_are_not_rigorous_error_bounds':True})
 assert abs(b['work_energy_error_Eh'])<1e-5 and abs(d['work_energy_error_Eh'])<1e-5
 # Independent tabulated reference has Ry energies: exactly .5Eh/Ry conversion.
 refs=json.loads((P/'data/independent_reference.json').read_text());checks=[]
 for row in refs['rows']:
  R=row['R_bohr'];gap=float(system.curves['Vu_Eh'](R)-system.curves['Vg_Eh'](R));checks.append({**row,'reference_gap_Eh':row['gap_Ry']/2,'computed_gap_Eh':gap,'difference_Eh':gap-row['gap_Ry']/2,'reference_Vg_Eh':(row['E0_Ry']-1)/2-row['gap_Ry']/4,'reference_Vu_Eh':(row['E0_Ry']-1)/2+row['gap_Ry']/4,'computed_Vg_Eh':float(system.curves['Vg_Eh'](R)),'computed_Vu_Eh':float(system.curves['Vu_Eh'](R))})
 save('independent_gap_checks.json',{'checks':checks,'source':refs['source'],'conversion':'1Ry = .5Eh explicitly; reference table approximations not an exact certified modern curve'})
 # Short backward compatibility only: old/new kernels on identical prepared state.
 import importlib.util
 spec=importlib.util.spec_from_file_location('old_core02',str(P/'data/core02_runtime_frozen.py'));old=importlib.util.module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)
 from h2_adapter import H2System
 from h3_adapter import H3System
 compat={}
 for name,model,init in [('H2',H2System(127),{'theta':.37,'relative_phase':.21}),('H3',H3System(128),{'energy_eV':.43,'center_R':5.7,'sigma_R':.72})]:
  w=QuantumWorld(model);w.prepare(init);oldstate=old.QuantumState(w.state.psi.copy(),model.identity)
  for _ in range(4):old.TimeEvolution.step(model,oldstate,.25)
  w.step(.25,count=4);error=float(np.max(abs(w.state.psi-oldstate.psi)));compat[name]={'four_step_state_error':error,'not_full_historical_experiment':True};assert error<1e-12
 save('Core02_compatibility.json',compat)
 B.checkpoint(OUT/'checkpoint');restored=QuantumWorld(system);restored.restore(OUT/'checkpoint');assert np.array_equal(restored.state.psi,B.state.psi) and restored.state.parameters==B.state.parameters
 probe=QuantumWorld(system);probe.prepare(initial);probe.apply(pulse);probe.step(.05,count=4)
 save('view_payload.json',{'model':system.metadata(),'model_id':system.identity,'allowed':system.allowed(),'R':system.R.tolist(),'Vg':system.hamiltonian.vg.tolist(),'Vu':system.hamiltonian.vu.tolist(),'mu':system.hamiltonian.mu.tolist(),'K':system.hamiltonian.K[0].tolist(),'dx':system.dx,'initial':initial,'pulse':B.state.parameters['pulse'],'actual_traces':{'A':a,'B':b,'C':c},'final_re':B.state.psi.real.tolist(),'final_im':B.state.psi.imag.tolist(),'browser_reference_re':probe.state.psi.real.tolist(),'browser_reference_im':probe.state.psi.imag.tolist()})
 subprocess.run(['node',str(P/'browser_qa.js'),str(OUT/'view_payload.json'),str(OUT/'browser_core_QA.json')],check=True,timeout=30)
 save('machine_QA.json',{'zero_field_decoupled':'PASS','real_field_transfer':'PASS','probability_balance':'PASS','fixed_R_analytic':'PASS','Core02_short_compatibility':'PASS','checkpoint':'PASS','physical_accuracy':'NOT_CERTIFIED','visual_QA':'NOT_RUN','classification':'PROGRAMMATIC_EXPERIMENT'})
 save('resources.json',{'wall_seconds':time.monotonic()-start,'cpu_user':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'cpu_system':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'cgroup':cgroup()})
if __name__=='__main__':
 try:main()
 except BaseException as exc:save('FAILURE.json',{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()});raise
