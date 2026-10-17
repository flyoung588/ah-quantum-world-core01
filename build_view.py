"""Assemble the standalone live view ONLY from actual cloud model/state output."""
import json
from pathlib import Path
p=Path(__file__).parent
payload=p/'results/view_payload.json'
if not payload.is_file():raise SystemExit('NOT_BUILT: actual cloud potential/state payload missing; no invented replacement')
data=payload.read_text();json.loads(data)
s=(p/'view-template.html').read_text().replace('__PAYLOAD__',data).replace('__ENGINE__',(p/'world-browser.js').read_text())
(p/'world.html').write_text(s)
