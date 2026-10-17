// Execute the actual browser physics core in Node ON THE SCIENCE RUNNER, not Mac.
const fs=require('fs'),Q=require('./world-browser.js'),p=JSON.parse(fs.readFileSync(process.argv[2]));
const w=new Q(p);let init=0;for(let i=0;i<w.re.length;i++)init=Math.max(init,Math.abs(w.re[i]-p.initial_re[i]),Math.abs(w.im[i]-p.initial_im[i]));
w.step(3);let error=0;for(let i=0;i<w.re.length;i++)error=Math.max(error,Math.abs(w.re[i]-p.reference_re[i]),Math.abs(w.im[i]-p.reference_im[i]));
const fork=w.fork();fork.apply('momentum_kick',{delta_k:-.75});fork.step(10);const comparison=w.compare(fork);
const qa={initial_max_error:init,three_step_Python_JS_max_error:error,pass:init<1e-12&&error<1e-11,fork_independent:w.t===1.5&&fork.t===6.5,comparison,export_has_state:fork.export_evidence().re.length===p.n*p.n,visual_browser_QA:'NOT_RUN; Node parity is not rendering QA'};
fs.writeFileSync(process.argv[3],JSON.stringify(qa,null,2));if(!qa.pass||!qa.fork_independent)process.exitCode=1;
