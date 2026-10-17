import base64,hashlib,json,pathlib,sys,zlib
root=pathlib.Path(sys.argv[1])
if root.is_dir():
 for f in sorted(root.rglob('*')):
  if not f.is_file():continue
  raw=f.read_bytes();s=base64.b64encode(zlib.compress(raw,9)).decode();parts=[s[i:i+7000] for i in range(0,len(s),7000)] or ['']
  for i,v in enumerate(parts):print('CORE02_FILE '+json.dumps({'name':str(f.relative_to(root)),'sha256':hashlib.sha256(raw).hexdigest(),'part':i,'parts':len(parts),'data':v}),flush=True)
