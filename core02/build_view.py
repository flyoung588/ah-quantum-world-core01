from pathlib import Path
import json,re
p=Path(__file__).parent
payload=p.parent/'results-core02/browser_models.json'
if not payload.exists():raise SystemExit('actual cloud output missing; no fake replacement')
data=json.loads(payload.read_text());data['H2']['validated_duration_au']=1200;data['H3']['validated_duration_au']=200
style=re.search(r'<style>(.*?)</style>',(p.parent/'view-template.html').read_text(),re.S).group(1)
html=(p/'view-template.html').read_text().replace('__STYLE__',style).replace('__PAYLOAD__',json.dumps(data,separators=(',',':'))).replace('__FFT__',(p.parent/'world-browser.js').read_text()).replace('__RUNTIME__',(p/'browser-runtime.js').read_text())
(p/'world.html').write_text(html)
