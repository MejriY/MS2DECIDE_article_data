from ms2decide.ClosestGNPS import _launch_GNPS_workflow_iterative
from ms2decide.AuthMail import AuthMail
from datetime import datetime
from pathlib import Path

input_dir = Path("../Generated/") / "Manufactured case/"
input_mgf = str((input_dir / "All GNPS.mgf").resolve())
input_quant = str((input_dir / "Quantification table.csv").resolve())
title = f"Manufactured case {datetime.now().isoformat()}"
auth = AuthMail.from_txt("../../../Auth GNPS.txt")
d = _launch_GNPS_workflow_iterative(auth, input_mgf, input_quant, title)
