from ms2decide.ClosestGNPS import _launch_GNPS_workflow_iterative
from ms2decide.ClosestGNPS import _upload_to_gnps
from ms2decide.AuthMail import AuthMail
from datetime import datetime
from pathlib import Path
import json

input_dir = Path("../Generated/") / "Manufactured case/"
input_mgf = str((input_dir / "All GNPS.mgf").resolve())
input_quant = str((input_dir / "Quantification table.csv").resolve())
auth = AuthMail.from_txt("../../../Auth GNPS.txt")
title = f"Manufactured case {datetime.now().isoformat()}"
(gnps_input_mgf, gnps_input_csv) = _upload_to_gnps(auth, input_mgf, input_quant, title)
d = _launch_GNPS_workflow_iterative(auth, gnps_input_mgf, gnps_input_csv, title)

# ids_by_first = {k:[v.split('=')[1] for kk,v in d.items()] for k,d in d_jobs.items()}
ids_by_first = {k:[v for kk,v in e.items()] for k,e in d.items()}
all_ids = [v for d in ids_by_first.values() for v in d]
with open(input_dir/"Gnps task ids.json", 'w') as f:
    f.write(json.dumps(all_ids))
