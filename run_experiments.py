"""One bounded cloud experiment, real synthetic-Hamiltonian and H3 TDSE tests."""
import json,os,time,traceback,resource,platform,subprocess
from pathlib import Path
import numpy as np
from world import *
OUT=Path('/out');OUT.mkdir(exist_ok=True)
def save(name,obj):
    p=OUT/name;tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False));tmp.replace(p)
def cg():
    p=Path('/sys/fs/cgroup');return {k:(p/k).read_text().strip() for k in ['memory.max','memory.current','memory.peak','memory.events','cpu.max','cpu.stat','memory.swap.max']}
def trajectory(e,energy,dt=.5,tmax=2000.,name=None):
    w=ExperimentControl(e);w.prepare(energy);rows=[w.observe()];snap=[]
    for i in range(int(tmax/(50*dt))):
        w.step(50,dt);rows.append(w.observe())
        if i%5==0:
            snap.append({'time_au':w.state.t,'density':(abs(w.state.psi[::4,::4])**2).tolist()})
    if name:
        save(name+'.json',{'initial':w.state.initial,'dt':dt,'n':e.n,'rows':rows,'snapshots_downsampled_only_for_display':snap})
        np.savez_compressed(OUT/(name+'-state.npz'),psi=w.state.psi,axis=e.axis)
    return w,rows
def main():
    start=time.monotonic();admit=cg()
    assert int(admit['memory.max'])<=1024**3 and admit['cpu.max']=='100000 100000' and admit['memory.swap.max']=='0'
    assert os.listdir('/sys/class/net')==['lo'] and len(os.sched_getaffinity(0))==1
    save('admission.json',{'cgroup':admit,'interfaces':os.listdir('/sys/class/net'),'affinity':list(os.sched_getaffinity(0))})
    save('environment.json',{'python':platform.python_version(),'numpy':np.__version__,'platform':platform.platform(),'sha':os.environ.get('GITHUB_SHA'),'mass_convention':'proton mass 1836.15267343 me, BO, not atomic mass'})
    e=PhysicsEngine(256)
    refs=[([1.757]*3,-.068262785352072364),([1.4,1000.,1001.4],-.17447419999999997),([1000.,1000.,2000.],-1.03502104795956e-17)]
    checks=[]
    for xyz,ref in refs:
        v=float(e.potential(*map(np.array,xyz)));checks.append({'xyz_bohr':xyz,'reference_Eh':ref,'actual_Eh':v,'error_Eh':v-ref})
    legacy=float(np.max(abs(e.V-e.potential(e.x,e.y,e.x+e.y,legacy=1))))
    saddle=float(e.potential(np.array(1.757),np.array(1.757),np.array(3.514)))
    save('pes_checks.json',{'upstream_points':checks,'legacy_refactoring_max_error_Eh':legacy,'collinear_saddle_Eh':saddle,'barrier_vs_upstream_min_eV':(saddle+.17447440045043025)*EV,'reference_barrier_eV_rounded':.425,'diatomic_energies_Eh':e.bond_energies.tolist()})
    assert max(abs(c['error_Eh']) for c in checks)<1e-10 and legacy<1e-12
    # Analytic 2D free Gaussian under the actual coupled kinetic operator.
    # Covariance evolves Sigma(t)=sigma^2 I+t^2 A A^T/(4 sigma^2), A=[[2,-1],[-1,2]]/m.
    a=e.axis;X,Y=e.x,e.y;sigma=.55;t=300.;p=np.array([.9,-.7]);center=np.array([6.,6.])
    psi=np.exp(-((X-6)**2+(Y-6)**2)/(4*sigma**2)+1j*(p[0]*X+p[1]*Y));psi/=np.linalg.norm(psi)*e.dx
    final=np.fft.ifft2(np.exp(-1j*e.K*t)*np.fft.fft2(psi));d=abs(final)**2*e.dx**2
    A=np.array([[2.,-1.],[-1.,2.]])/MASS;expected_center=center+t*A@p;expected_cov=sigma**2*np.eye(2)+t*t*A@A.T/(4*sigma**2)
    mean=np.array([(d*X).sum(),(d*Y).sum()]);cov=np.array([[(d*(X-mean[0])**2).sum(),(d*(X-mean[0])*(Y-mean[1])).sum()],[(d*(X-mean[0])*(Y-mean[1])).sum(),(d*(Y-mean[1])**2).sum()]])
    save('analytic_free_packet.json',{'actual_center':mean.tolist(),'analytic_center':expected_center.tolist(),'actual_cov':cov.tolist(),'analytic_cov':expected_cov.tolist(),'center_max_error':float(abs(mean-expected_center).max()),'cov_max_error':float(abs(cov-expected_cov).max()),'norm':float(d.sum())})
    assert abs(mean-expected_center).max()<1e-8 and abs(cov-expected_cov).max()<1e-8
    # Closed H3 propagation: no absorber. Check norm and energy, not CAP conservation.
    closed=PhysicsEngine(256,cap=False);cw,cr=trajectory(closed,.6,tmax=400.)
    closedqa={'norm_max_error':max(abs(v['norm']-1) for v in cr),'energy_max_drift_Eh':max(abs(v['energy_Eh']-cr[0]['energy_Eh']) for v in cr),'max_edge_probability':max(v['edge_probability'] for v in cr),'rows':cr}
    save('closed_system.json',closedqa)
    assert closedqa['norm_max_error']<1e-9
    assert closedqa['energy_max_drift_Eh']<1e-6
    # Three new allowed preparations, no static H3 data imported.
    worlds={};runs={}
    for en in [.2,.6,1.]:
        w,rows=trajectory(e,en,name='experiment-'+str(en));worlds[en]=w;runs[en]=rows
    # Real feedback: compare lower/higher energy, then apply a bounded kick to a fork.
    initial=ExperimentControl(e);initial.prepare(.2);branch=initial.fork();before=initial.observe();branch.apply('momentum_kick',delta_k=-.75)
    initial.step(800);branch.step(800)
    save('ai_loop.json',{'request':'change incident momentum after seeing the lower/higher-energy comparison','energy_comparison':worlds[.2].compare(worlds[1.]),'fork_comparison':initial.compare(branch),'operation':branch.state.controls,'fork_preserved_at_creation':before,'new_observations_not_static_lookup':True})
    # Independent initial-state regeneration on refined grids; compare regional populations, not mismatched wavefunctions.
    fine_dt,rd=trajectory(e,.6,dt=.25,name='refine-dt')
    fine_grid,rg=trajectory(PhysicsEngine(512),.6,name='refine-grid')
    base=runs[.6][-1];dtlast=rd[-1];glast=rg[-1]
    conv={k:{'baseline':base[k],'dt_halved':dtlast[k],'grid_doubled':glast[k],'dt_difference':abs(base[k]-dtlast[k]),'grid_difference':abs(base[k]-glast[k])} for k in ['product_region','reactant_region','absorbed','energy_Eh']}
    conv['criteria']={'dt_product_abs':1e-3,'grid_product_abs':1e-2,'meaning':'limited regional-population convergence, not state-to-state reaction S-matrix'}
    conv['within_frozen_regional_tolerances']=bool(conv['product_region']['dt_difference']<1e-3 and conv['product_region']['grid_difference']<1e-2)
    save('convergence.json',conv)
    capqa={str(en):{'max_probability_balance_error':max(abs(v['probability_balance']-1) for v in rows),'max_edge_probability':max(v['edge_probability'] for v in rows),'final':rows[-1]} for en,rows in runs.items()}
    save('cap_accounting.json',capqa)
    assert max(v['max_probability_balance_error'] for v in capqa.values())<1e-9
    # Real Python state exported to browser engine; view computes subsequent states, no prerecorded animation.
    view=ExperimentControl(e);view.prepare(.6);initialpsi=view.state.psi.copy();view.step(3)
    save('view_payload.json',{'n':e.n,'dx':e.dx,'axis':e.axis.tolist(),'V':e.V.ravel().tolist(),'W':e.W.ravel().tolist(),'K':e.K.ravel().tolist(),'bond':e.bond.tolist(),'mass':MASS,'ev':EV,'dt':.5,'reference_re':view.state.psi.real.ravel().tolist(),'reference_im':view.state.psi.imag.ravel().tolist(),'initial_re':initialpsi.real.ravel().tolist(),'initial_im':initialpsi.imag.ravel().tolist()})
    subprocess.run(['node','browser_qa.js',str(OUT/'view_payload.json'),str(OUT/'browser_engine_qa.json')],check=True,timeout=40)
    save('machine_qa.json',{'PES_reference':'PASS','analytic_free_wavepacket':'PASS','closed_norm':'PASS','closed_energy_drift_Eh':closedqa['energy_max_drift_Eh'],'CAP_probability_accounting':'PASS','regional_convergence':conv['within_frozen_regional_tolerances'],'physics_claim':'approximate collinear single-surface LSTH TDSE, no new physical discovery'})
    save('resources.json',{'wall_seconds':time.monotonic()-start,'cpu_user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'cpu_system_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'ru_maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'cgroup':cg()})
if __name__=='__main__':
    try:main()
    except BaseException as exc:
        save('FAILURE.json',{'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()});raise
