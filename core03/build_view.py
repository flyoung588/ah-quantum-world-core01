"""Build only after actual physical output exists; no substitute surfaces/dipoles."""
import json,re
from pathlib import Path
p=Path(__file__).parent;r=p.parent;new=r/'results-core03/view_payload.json'
if not new.exists():raise SystemExit('NOT_BUILT: actual H2+ electronic and time-evolution payload missing')
models=json.loads((r/'results-core02/browser_models.json').read_text());newdata=json.loads(new.read_text())
newdata.update({'metadata':newdata['model'],'model_id':newdata['model_id'],'dv':newdata['dx'],'axes':[newdata['R']],'allowed':newdata['allowed'],'dt_au':.05,'validated_duration_au':240})
models['H2plus']=newdata;models['H2']['validated_duration_au']=1200;models['H3']['validated_duration_au']=200
html=(r/'core02/view-template.html').read_text().replace('Core-02','Core-03')
html=html.replace('</select>','<option value="H2plus">H₂⁺ Coupled Electronic Channels</option></select>')
html=html.replace('makers={H2:h2Adapter,H3:h3Adapter}','makers={H2:h2Adapter,H3:h3Adapter,H2plus:h2plusAdapter}')
html=html.replace('<details>','<label>允许控制 JSON <textarea id="control" rows="5" style="width:95%"></textarea></label><details>',1)
html=html.replace("$('input').value=JSON.stringify(world.p.initial,null,2);","$('input').value=JSON.stringify(world.p.initial,null,2);$('control').value=JSON.stringify(world.adapter.kick,null,2);")
html=html.replace('world.apply(world.adapter.kick)',"world.apply(JSON.parse($('control').value))")
html=html.replace('world.step(.25,4)','world.step(world.p.dt_au||.25,4)').replace('world.step(.25,2)','world.step(world.p.dt_au||.25,2)')
html=html.replace('模型允许的相位操作','施加模型允许控制').replace('height="400"','height="600"')
style=re.search(r'<style>(.*?)</style>',(r/'view-template.html').read_text(),re.S).group(1)
html=html.replace('__STYLE__',style).replace('__PAYLOAD__',json.dumps(models,separators=(',',':'))).replace('__FFT__',(r/'world-browser.js').read_text()).replace('__RUNTIME__',(r/'core02/browser-runtime.js').read_text()+'\n'+(p/'browser-adapter.js').read_text())
(p/'world.html').write_text(html)
