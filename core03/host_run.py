import json,os,pathlib,subprocess,sys
out=pathlib.Path(os.environ['RUNNER_TEMP'])/'core03-results';out.mkdir(exist_ok=True)
first=min(os.sched_getaffinity(0))
cmd=['docker','create','--network','none','--memory','1g','--memory-swap','1g','--cpus','1','--cpuset-cpus',str(first),'--pids-limit','64','--cap-drop','ALL','--security-opt','no-new-privileges','--read-only','--tmpfs','/tmp:rw,noexec,nosuid,size=64m','--user',f'{os.getuid()}:{os.getgid()}','--mount',f'type=bind,src={out},dst=/out']
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:cmd+=['-e',k+'=1']
cmd+=['-e','GITHUB_SHA='+os.environ['GITHUB_SHA'],'core03:once'];cid=subprocess.check_output(cmd,text=True).strip();code=1
try:
 info=json.loads(subprocess.check_output(['docker','inspect',cid]))[0];h=info['HostConfig']
 assert h['Memory']==1024**3 and h['MemorySwap']==1024**3 and h['NetworkMode']=='none' and h['NanoCpus']==1000000000 and not h['Privileged']
 (out/'host_container.json').write_text(json.dumps({'HostConfig':h,'image':info['Image'],'uid':info['Config']['User']},indent=2))
 code=subprocess.run(['docker','start','-a',cid],timeout=730).returncode
except BaseException as exc:
 (out/'host_failure.json').write_text(json.dumps({'type':type(exc).__name__,'message':str(exc)}))
finally:
 subprocess.run(['docker','kill',cid],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=3)
 subprocess.run(['docker','rm',cid],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=3)
sys.exit(code)
