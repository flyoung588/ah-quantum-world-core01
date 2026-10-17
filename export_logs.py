"""Export exact bytes through free ordinary logs; no artifact/cache storage."""
import base64,hashlib,json,pathlib,sys,zlib
for p in sorted(pathlib.Path(sys.argv[1]).iterdir()):
    if not p.is_file() or p.name.startswith('.'):continue
    raw=p.read_bytes();s=base64.b64encode(zlib.compress(raw,9)).decode();parts=[s[i:i+7000] for i in range(0,len(s),7000)] or ['']
    for i,v in enumerate(parts):print('CORE_FILE '+json.dumps({'name':p.name,'sha256':hashlib.sha256(raw).hexdigest(),'part':i,'parts':len(parts),'data':v}),flush=True)
