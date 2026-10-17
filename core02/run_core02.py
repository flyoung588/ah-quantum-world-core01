"""One bounded cross-system experiment via identical runtime flow. No Core01 rerun."""
import hashlib,json,os,platform,resource,subprocess,time,traceback,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,scipy
from runtime import QuantumWorld,clean
from h2_adapter import H2System,CM
from h3_adapter import H3System
OUT=Path('/out');OUT.mkdir(exist_ok=True)
ROOT=Path(__file__).parent
def save(name,obj):
    p=OUT/name;tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(clean(obj),indent=2,allow_nan=False));tmp.replace(p)
def cg():return {k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in ['memory.max','memory.peak','memory.current','memory.events','memory.swap.max','cpu.max','cpu.stat']}
def run(system,request,duration=None,dt=.25,record_density=False):
    # Identical control flow for any PhysicalSystem. No identity switch here.
    w=QuantumWorld(system);allowed=w.allowed();w.prepare(request['initial']);initial=w.fork();rows=[w.observe()];density=[]
    duration=request['duration_au'] if duration is None else duration
    for i in range(round(duration/(10*dt))):
        w.step(dt,count=10);rows.append(w.observe())
        if record_density and i%6==0:density.append({'time_au':w.state.t,'density':(abs(w.state.psi)**2).ravel().tolist()})
    branch=initial.fork();branch.apply(request['control']);untouched=initial.fork()
    branch.evolve(request['branch_duration_au'],dt);untouched.evolve(request['branch_duration_au'],dt)
    comparison=untouched.compare(branch)
    return w,rows,{'allowed':allowed,'input':request,'observations':rows,'density':density,'branch_comparison':comparison,'control_work':branch.state.controls,'selection_status':'PROGRAMMATIC_EXPERIMENT'}
def browser_payload(world):
    s=world.system;h=s.hamiltonian;meta=s.metadata()
    common={'metadata':meta,'model_id':s.identity,'shape':list(h.V.shape),'dv':s.dv,'V':h.V.ravel().tolist(),'K':h.K.ravel().tolist(),'W':h.W.ravel().tolist(),'axes':[x.tolist() for x in s.axes],'allowed':s.allowed()}
    return common
def main():
    start=time.monotonic();admit=cg()
    assert admit['memory.max']=='1073741824' and admit['memory.swap.max']=='0' and admit['cpu.max']=='100000 100000'
    assert os.listdir('/sys/class/net')==['lo'] and len(os.sched_getaffinity(0))==1
    save('admission.json',{'cgroup':admit,'network':os.listdir('/sys/class/net'),'affinity':list(os.sched_getaffinity(0))})
    save('environment.json',{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sha':os.environ.get('GITHUB_SHA'),'threads':{k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}})
    requests=json.loads((ROOT/'requests.json').read_text());save('actual_requests.json',requests)
    # Prove H2 neither imports nor opens H3/LSTH: isolated subprocess with library hidden.
    check=subprocess.run([sys.executable,'-c',"import sys;from h2_adapter import H2System;h=H2System(127);assert 'world' not in sys.modules and 'h3_adapter' not in sys.modules;print(h.metadata())"],capture_output=True,text=True,timeout=15)
    save('h2_dependency_check.json',{'exit_code':check.returncode,'stdout':check.stdout,'stderr':check.stderr,'meaning':'H2 module and potential instantiate without importing H3 or loading ctypes library'});assert check.returncode==0
    h2=H2System();h3=H3System()
    # Frozen Core01 actual reference: initial state, three-step state, and operator hashes.
    ref=json.loads((ROOT/'data/core01_3step.json').read_text())
    regression=QuantumWorld(h3);regression.prepare({'energy_eV':.6,'center_R':6.,'sigma_R':.6})
    initial_error=float(np.max(abs(regression.state.psi-(np.array(ref['initial_re'])+1j*np.array(ref['initial_im'])).reshape(256,256))))
    regression.step(.5,count=3)
    state_error=float(np.max(abs(regression.state.psi-(np.array(ref['reference_re'])+1j*np.array(ref['reference_im'])).reshape(256,256))))
    hashes={k:hashlib.sha256(json.dumps(getattr(h3.engine,k).ravel().tolist(),separators=(',',':')).encode()).hexdigest() for k in ['V','K','W']}
    save('H3_frozen_regression.json',{'initial_max_error':initial_error,'three_step_state_error':state_error,'operator_hashes':hashes,'frozen_operator_hashes':ref['operator_hashes'],'full_experiment_rerun':False,'mixed_kinetic':'(kx²+ky²-kx*ky)/m unchanged'})
    assert initial_error<1e-12 and state_error<1e-11 and hashes==ref['operator_hashes']
    results={};worlds={};browser={}
    # Same run() implementation, different declared adapters and run-time JSON.
    for name,system in [('H2',h2),('H3',h3)]:
        w,rows,result=run(system,requests[name],record_density=name=='H2');worlds[name]=w;results[name]=result
        save(name+'_experiment.json',result)
        w.error_evidence={'local_probability_balance_max_error':max(abs(x['probability_balance']-1) for x in rows),'local_energy_max_drift_Eh':max(abs(x['energy_Eh']-rows[0]['energy_Eh']) for x in rows),'convergence_report':name+'_convergence_and_physics.json' if name=='H2' else 'H3_short_convergence.json','meaning':'local numerical diagnostics, not independent physical certification'}
        cp=OUT/(name+'_checkpoint');w.checkpoint(cp)
        restored=QuantumWorld(system);restored.restore(cp);assert np.array_equal(restored.state.psi,w.state.psi)
        w.export_evidence(OUT/(name+'_evidence'))
        save(name+'_QA.json',{'checkpoint_exact':True,'norm_error':max(abs(x['norm']-1) for x in rows),'probability_balance_error':max(abs(x['probability_balance']-1) for x in rows),'energy_drift_Eh':max(abs(x['energy_Eh']-rows[0]['energy_Eh']) for x in rows),'branch_same_time':result['branch_comparison']['same_time'],'branch_state_L2':result['branch_comparison']['state_L2']})
        assert max(abs(x['probability_balance']-1) for x in rows)<1e-9
        assert max(abs(x['energy_Eh']-rows[0]['energy_Eh']) for x in rows)<1e-6
        assert result['branch_comparison']['state_L2']>1e-6 and result['branch_comparison']['same_time']
        probe=QuantumWorld(system);probe.prepare(requests[name]['initial']);payload=browser_payload(probe)
        if name=='H2':payload.update({'transform':'sine','basis0':h2.eigenvectors[:,0].tolist(),'basis1':h2.eigenvectors[:,1].tolist(),'initial':requests[name]['initial']})
        else:payload.update({'transform':'fft2','bond':h3.engine.bond.tolist(),'mass':1836.15267343,'initial':requests[name]['initial']})
        probe.step(.25,count=4);payload['reference_re']=probe.state.psi.real.ravel().tolist();payload['reference_im']=probe.state.psi.imag.ravel().tolist();browser[name]=payload
    rejected=False
    try:worlds['H2'].compare(worlds['H3'])
    except ValueError:rejected=True
    assert rejected
    save('cross_system_contract.json',{'same_runtime_class':True,'same_main_control_flow':'run(system,request) for both','cross_model_subtraction_rejected':rejected,'models':{k:w.system.metadata() for k,w in worlds.items()},'not_universal':'grid states and diagonal-potential spectral split operators only'})
    # H2 numerical errors and potential-interpolation errors are deliberately separated.
    fine,fr,_=run(h2,requests['H2'],dt=.125)
    dense,gr,_=run(H2System(511),requests['H2'])
    extended,er,_=run(H2System(511,lower=.1,upper=10.),requests['H2'])
    alternative=H2System(interpolation='pchip')
    gap=float((h2.levels[1]-h2.levels[0])*CM);gap_grid=float((dense.system.levels[1]-dense.system.levels[0])*CM)
    gap_alt=float((alternative.levels[1]-alternative.levels[0])*CM)
    last=results['H2']['observations'][-1]
    convergence={'time_step_meanR_difference':abs(last['mean_R']-fr[-1]['mean_R']),'grid_meanR_difference':abs(last['mean_R']-gr[-1]['mean_R']),'extended_box_meanR_difference':abs(last['mean_R']-er[-1]['mean_R']),'gap_cm_inverse':gap,'grid_gap_difference_cm_inverse':gap_grid-gap,'pchip_gap_cm_inverse':gap_alt,'interpolation_sensitivity_cm_inverse':gap_alt-gap,'interpolation_sensitivity_is_NOT_rigorous_error_bound':True,'independent_BO_reference_cm_inverse':4163.40350,'independent_experimental_reference_cm_inverse':4161.16632,'BO_discrepancy_cm_inverse':gap-4163.40350,'experiment_discrepancy_cm_inverse':gap-4161.16632,'reference':'Dickenson etal PRL110193601(2013) table2 https://arxiv.org/html/1304.6390','missing_physics':'adiabatic, nonadiabatic, relativistic, QED terms absent; experiment mismatch NOT pure numerical error'}
    # Actual time-dependent bond-length signal vs exact evolution of same discrete H.
    probe=QuantumWorld(h2);probe.prepare(requests['H2']['initial']);psi0=probe.state.psi;exact=h2.eigenvectors@((h2.eigenvectors.T@psi0)*np.exp(-1j*h2.levels*1200))
    convergence['state_L2_vs_exact_same_discrete_H']=float(np.linalg.norm(worlds['H2'].state.psi-exact)*np.sqrt(h2.dx))
    convergence['same_H_check_not_independent_physics']=True
    save('H2_convergence_and_physics.json',convergence)
    # Short H3 dt/grid checks only; no48fs/full3-energy reruns.
    dt3,dr3,_=run(h3,requests['H3'],dt=.125)
    coarse,cr3,_=run(H3System(128),requests['H3'])
    save('H3_short_convergence.json',{'t_au':200,'dt_mean_y_difference':abs(results['H3']['observations'][-1]['mean_y']-dr3[-1]['mean_y']),'grid128_vs256_mean_y_difference':abs(results['H3']['observations'][-1]['mean_y']-cr3[-1]['mean_y']),'not_asymptotic_reaction_probability':True})
    # CAP diagnostic in H2: same run-time control, nonzero exterior tail attenuated.
    absorbed=QuantumWorld(H2System(cap=True));absorbed.prepare({'kind':'radial_gaussian','center':6.8,'sigma':.25,'k':.5});absorbed.evolve(40.,.25)
    save('H2_CAP_accounting.json',absorbed.observe())
    assert absorbed.observe()['absorbed']>1e-4 and abs(absorbed.observe()['probability_balance']-1)<1e-9
    save('browser_models.json',browser)
    subprocess.run(['node',str(ROOT/'browser_core_qa.js'),str(OUT/'browser_models.json'),str(OUT/'browser_core_QA.json')],check=True,timeout=30)
    save('machine_QA.json',{'cross_system_runtime':'PASS','H3_frozen_regression':'PASS','H2_independent_dependency':'PASS','checkpoint_restore':'PASS','cross_system_subtraction_rejected':'PASS','actual_data_not_static_lookup':True,'H2_precision_spectroscopy':'NOT_QUALIFIED','selection':'PROGRAMMATIC_EXPERIMENT, not AUTONOMOUS_DISCOVERY','visual_browser_QA':'NOT_RUN'})
    save('resources.json',{'wall_seconds':time.monotonic()-start,'cpu_user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'cpu_system_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'cgroup':cg()})
if __name__=='__main__':
 try:main()
 except BaseException as exc:
  save('FAILURE.json',{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()});raise
