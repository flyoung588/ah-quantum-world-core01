"""H2 J=0 reduced radial wavefunction u(R); independent of ALL H3/LSTH code."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.fft import dst,idst
from scipy.interpolate import CubicSpline,PchipInterpolator
from runtime import PhysicalSystem,Hamiltonian
MASS=1836.15267343;CM=219474.6313632

class H2System(PhysicalSystem):
    max_dt=.5
    def __init__(self,n=255,lower=.2,upper=8.,interpolation='cubic',cap=False):
        if n not in [127,255,511] or lower<.1 or upper>20 or lower>=upper:raise ValueError('radial model domain')
        self.n=n;self.lower=lower;self.upper=upper;self.interpolation=interpolation;self.cap=cap
        self.dx=(upper-lower)/(n+1);self.R=lower+self.dx*np.arange(1,n+1);self.axes=[self.R];self.dv=self.dx
        datafile=Path(__file__).parent/'data/h2_bo.json';self.data_sha=hashlib.sha256(datafile.read_bytes()).hexdigest();table=json.loads(datafile.read_text());r,e=np.array(table['rows']).T
        self.potential=CubicSpline(r,e,extrapolate=False) if interpolation=='cubic' else PchipInterpolator(r,e,extrapolate=False)
        V=self.potential(self.R);K=(np.pi*np.arange(1,n+1)/(upper-lower))**2/MASS
        W=.02*np.maximum((self.R-6.)/2.,0)**4 if cap else np.zeros(n)
        forward=lambda p:dst(p,type=1,norm='ortho');inverse=lambda p:idst(p,type=1,norm='ortho')
        self.hamiltonian=Hamiltonian(V,K,forward,inverse,W,'symmetric split operator + DST-I sine-DVR; T=-dRR/m, mu=mp/2')
        S=dst(np.eye(n),type=1,norm='ortho',axis=0);matrix=(S*K[None,:])@S.T+np.diag(V)
        self.levels,self.eigenvectors=np.linalg.eigh(matrix)
        # Deterministic eigenvector signs make prepared phases reproducible across grids.
        for j in range(n):
            i=np.argmax(abs(self.eigenvectors[:,j]));self.eigenvectors[:,j]*=np.sign(self.eigenvectors[i,j])
    def metadata(self):return {'system':'H2 J=0 nuclear vibration','particles':{'protons':2,'electrons':2,'nuclear_charge':1,'proton_mass_me':MASS,'reduced_mass_me':MASS/2,'electronic_spin':'singlet ground X1Sigma_g+; nuclear para-H2 J=0'},'coordinates':'R>0, reduced radial u(R); probability |u|² dR, NOT R²|u|²dR','units':'bohr,Eh,atomic time','energy_zero':'two separated neutral H atoms=-1Eh; TOTAL BO energy convention','source':'Pachucki PRA82 032509(2010) table4; 82 published numerical pairs; '+self.data_sha,'interpolation':self.interpolation,'boundary':{'type':'Dirichlet u(lower)=u(upper)=0; no periodic wrap','lower':self.lower,'upper':self.upper,'CAP':self.cap},'n':self.n,'approximation':'single BO surface J=0; no adiabatic/nonadiabatic/relativistic/QED corrections','scope':'table-interpolated model dynamics, NOT precision spectroscopy; interpolation/box error separately measured'}
    def allowed(self):return {'initial':{'theta':[0.,.8],'relative_phase':[-float(np.pi),float(np.pi)],'meaning':'cos(theta)|v0>+exp(i phase)sin(theta)|v1> of this discrete Hamiltonian','alternative_radial_gaussian':{'center':[1.,7.],'sigma':[.2,.5],'k':[-1.,1.]}},'controls':{'radial_phase_imprint':{'q_bohr_inverse':[-1.,1.],'meaning':'ideal unitary exp(i q R), NOT literal laser control; work recorded'}},'observables':['mean_R','variance_R','vibrational_populations','edge_probability']}
    def initial_wavefunction(self,initial):
        if initial.get('kind')=='radial_gaussian':
            c=float(initial['center']);sigma=float(initial['sigma']);k=float(initial['k'])
            if not(np.isfinite([c,sigma,k]).all() and 1<=c<=7 and .2<=sigma<=.5 and abs(k)<=1):raise ValueError('Gaussian bounds')
            return np.exp(-(self.R-c)**2/(4*sigma*sigma)+1j*k*self.R)*np.sin(np.pi*(self.R-self.lower)/(self.upper-self.lower))
        theta=float(initial['theta']);phase=float(initial['relative_phase'])
        if not(np.isfinite(theta) and np.isfinite(phase) and 0<=theta<=.8 and abs(phase)<=np.pi):raise ValueError('H2 preparation bounds')
        return (np.cos(theta)*self.eigenvectors[:,0]+np.exp(1j*phase)*np.sin(theta)*self.eigenvectors[:,1])/np.sqrt(self.dx)
    def apply_control(self,psi,control):
        if control.get('operation')!='radial_phase_imprint':raise ValueError('unsupported H2 control')
        q=float(control['q_bohr_inverse'])
        if not np.isfinite(q) or abs(q)>1:raise ValueError('phase bound')
        return psi*np.exp(1j*q*self.R)
    def observables(self,state):
        d=abs(state.psi)**2*self.dx;norm=d.sum();mean=float((d*self.R).sum()/norm)
        pop=abs(self.eigenvectors[:,:6].T@state.psi*np.sqrt(self.dx))**2
        return {'mean_R':mean,'variance_R':float((d*(self.R-mean)**2).sum()/norm),'vibrational_populations':pop.tolist(),'edge_probability':float(d[(self.R<.5)|(self.R>7.)].sum())}
