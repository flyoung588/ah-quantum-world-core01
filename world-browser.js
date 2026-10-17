/* Same spectral split-operator as world.py. No network, no prerecorded frames. */
'use strict';
class QuantumWorld {
  constructor(p){this.p=p;this.n=p.n;this.re=new Float64Array(p.n*p.n);this.im=new Float64Array(p.n*p.n);this.prepare();}
  prepare(energy=.6,center=6,sigma=.6){
    if(!(energy>=.05&&energy<=1.2&&center>=4.5&&center<=8&&sigma>=.45&&sigma<=1))throw Error('outside allowed preparation domain');
    const p=this.p,n=this.n,k=-Math.sqrt(2*(2*p.mass/3)*energy/p.ev);let norm=0;
    for(let i=0;i<n;i++)for(let j=0;j<n;j++){const q=i*n+j,R=p.axis[j]+p.axis[i]/2,a=p.bond[i]*Math.exp(-((R-center)**2)/(4*sigma*sigma));this.re[q]=a*Math.cos(k*R);this.im[q]=a*Math.sin(k*R);norm+=a*a*p.dx*p.dx;}
    const f=1/Math.sqrt(norm);for(let q=0;q<n*n;q++){this.re[q]*=f;this.im[q]*=f;}
    this.t=0;this.absorbed=0;this.initial={energy_eV:energy,center_R:center,sigma_R:sigma};this.controls=[];return this.observe();
  }
  static fft(r,im,inverse){
    const n=r.length;for(let i=1,j=0;i<n;i++){let b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;if(i<j){[r[i],r[j]]=[r[j],r[i]];[im[i],im[j]]=[im[j],im[i]];}}
    for(let len=2;len<=n;len<<=1){const a=(inverse?2:-2)*Math.PI/len,wr0=Math.cos(a),wi0=Math.sin(a);for(let i=0;i<n;i+=len){let wr=1,wi=0;for(let j=0;j<len/2;j++){const u=i+j,v=u+len/2,vr=r[v]*wr-im[v]*wi,vi=r[v]*wi+im[v]*wr;r[v]=r[u]-vr;im[v]=im[u]-vi;r[u]+=vr;im[u]+=vi;const t=wr*wr0-wi*wi0;wi=wr*wi0+wi*wr0;wr=t;}}}
    if(inverse)for(let i=0;i<n;i++){r[i]/=n;im[i]/=n;}
  }
  fft2(inv){const n=this.n,r=new Float64Array(n),im=new Float64Array(n);for(let i=0;i<n;i++){r.set(this.re.subarray(i*n,(i+1)*n));im.set(this.im.subarray(i*n,(i+1)*n));QuantumWorld.fft(r,im,inv);this.re.set(r,i*n);this.im.set(im,i*n);}for(let j=0;j<n;j++){for(let i=0;i<n;i++){r[i]=this.re[i*n+j];im[i]=this.im[i*n+j];}QuantumWorld.fft(r,im,inv);for(let i=0;i<n;i++){this.re[i*n+j]=r[i];this.im[i*n+j]=im[i];}}}
  phase(c,s){for(let q=0;q<this.re.length;q++){const r=this.re[q],im=this.im[q];this.re[q]=r*c[q]-im*s[q];this.im[q]=r*s[q]+im*c[q];}}
  step(count=1,dt=.5){
    if(!(Number.isInteger(count)&&count>=1&&count<=500&&dt>0&&dt<=.5))throw Error('bounded step');
    const p=this.p,m=this.re.length,c=new Float64Array(m),s=new Float64Array(m),kc=new Float64Array(m),ks=new Float64Array(m);
    for(let q=0;q<m;q++){const a=-p.V[q]*dt/2,decay=Math.exp(-p.W[q]*dt/2);c[q]=decay*Math.cos(a);s[q]=decay*Math.sin(a);kc[q]=Math.cos(p.K[q]*dt);ks[q]=-Math.sin(p.K[q]*dt);}
    const loss=new Float64Array(m);for(let q=0;q<m;q++)loss[q]=-Math.expm1(-p.W[q]*dt);
    for(let k=0;k<count;k++){let absorbed=0;for(let q=0;q<m;q++)absorbed+=(this.re[q]**2+this.im[q]**2)*loss[q]*p.dx**2;this.phase(c,s);this.fft2(false);this.phase(kc,ks);this.fft2(true);for(let q=0;q<m;q++)absorbed+=(this.re[q]**2+this.im[q]**2)*loss[q]*p.dx**2;this.phase(c,s);this.absorbed+=absorbed;this.t+=dt;}
    return this.observe();
  }
  norm(){let n=0;for(let q=0;q<this.re.length;q++)n+=this.re[q]**2+this.im[q]**2;return n*this.p.dx**2;}
  observe(){let norm=0,mx=0,my=0,product=0,reactant=0,interaction=0;const n=this.n,p=this.p;for(let i=0;i<n;i++)for(let j=0;j<n;j++){const q=i*n+j,d=(this.re[q]**2+this.im[q]**2)*p.dx**2,x=p.axis[i],y=p.axis[j];norm+=d;mx+=x*d;my+=y*d;if(x>3&&y<3)product+=d;if(x<3&&y>3)reactant+=d;if(x<3&&y<3)interaction+=d;}return{time_au:this.t,time_fs:this.t*.024188843265857,norm,absorbed:this.absorbed,probability_balance:norm+this.absorbed,mean_x:mx/norm,mean_y:my/norm,product_region:product,reactant_region:reactant,interaction_region:interaction};}
  apply(operation,{delta_k}={}){if(operation!=='momentum_kick')throw Error('not implemented');if(!Number.isFinite(delta_k)||Math.abs(delta_k)>1)throw Error('bounded kick');const n=this.n,p=this.p;for(let i=0;i<n;i++)for(let j=0;j<n;j++){const q=i*n+j,a=delta_k*(p.axis[j]+p.axis[i]/2),r=this.re[q],im=this.im[q];this.re[q]=r*Math.cos(a)-im*Math.sin(a);this.im[q]=r*Math.sin(a)+im*Math.cos(a);}this.controls.push({operation,delta_k,time_au:this.t});return this.observe();}
  fork(){const b=Object.create(QuantumWorld.prototype);b.p=this.p;b.n=this.n;b.re=this.re.slice();b.im=this.im.slice();b.t=this.t;b.absorbed=this.absorbed;b.initial={...this.initial};b.controls=JSON.parse(JSON.stringify(this.controls));return b;}
  compare(other){return{first:this.observe(),second:other.observe(),product_difference:other.observe().product_region-this.observe().product_region};}
  export_evidence(){return{model:'LSTH collinear TDSE; finite grid',initial:this.initial,controls:this.controls,observation:this.observe(),n:this.n,dx:this.p.dx,re:Array.from(this.re),im:Array.from(this.im),ceiling:'regional populations, not asymptotic reaction rates; browser state unqualified beyond cloud convergence checks'};}
}
if(typeof module!=='undefined')module.exports=QuantumWorld;
