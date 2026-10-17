"""Physical H3 semantics moved out of runtime; original immutable model reused."""
import numpy as np
from runtime import PhysicalSystem,Hamiltonian
from world import PhysicsEngine,ExperimentControl,ObservableEngine,MASS

class H3System(PhysicalSystem):
    max_dt=.5
    def __init__(self,n=256,cap=True):
        self.engine=PhysicsEngine(n,cap=cap);e=self.engine;self.n=n;self.cap=cap
        self.axes=[e.axis,e.axis];self.dv=e.dx**2
        self.hamiltonian=Hamiltonian(e.V,e.K,lambda p:np.fft.fftn(p,norm='ortho'),lambda p:np.fft.ifftn(p,norm='ortho'),e.W,'symmetric split operator + 2D FFT; mixed kinetic kx²+ky²-kxky /m')
    def metadata(self):return {'system':'neutral H3 collinear','particles':{'protons':3,'electrons':3,'nuclear_charge':1,'mass_me':MASS,'electronic_spin':'doublet ground surface'},'coordinates':'x=rAB,y=rBC, ordered nuclei','units':'bohr,Eh,atomic time','energy_zero':'three separated neutral H atoms in LSTH convention','source':'QML LSTH commit73f2ef3d322b253a4ef39e08216261bd974180f6; MIT','boundary':'periodic spectral box .4..14.4bohr; short-range repulsion; '+('large-separation CAP' if self.cap else 'no CAP'),'n':self.n,'approximation':'collinear single adiabatic BO surface, no rotational/nonadiabatic dynamics','scope':'short regression only; no new asymptotic reaction probability'}
    def allowed(self):return {'initial':{'energy_eV':[.05,1.2],'center_R':[4.5,8.],'sigma_R':[.45,1.]},'controls':{'momentum_kick':{'delta_k':[-1.,1.],'coordinate':'R=y+x/2','implementation':'ideal instantaneous phase imprint; work recorded'}},'observables':['mean_x','mean_y','reactant_region','product_region','interaction_region','edge_probability']}
    def initial_wavefunction(self,initial):
        old=ExperimentControl(self.engine);old.prepare(**initial);return old.state.psi
    def apply_control(self,psi,control):
        if control.get('operation')!='momentum_kick':raise ValueError('unsupported H3 control')
        k=float(control['delta_k'])
        if not np.isfinite(k) or abs(k)>1:raise ValueError('invalid kick')
        return psi*np.exp(1j*k*(self.engine.y+self.engine.x/2))
    def observables(self,state):
        o=ObservableEngine.measure(self.engine,state)
        return {k:v for k,v in o.items() if k in self.allowed()['observables']}
