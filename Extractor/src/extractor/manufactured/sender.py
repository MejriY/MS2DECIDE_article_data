from ms2decide.ClosestGNPS import _invoke_workflow
from ms2decide.ClosestGNPS import _get_networking_parameters
from ms2decide.ClosestGNPS import _upload_to_gnps
from ms2decide.AuthMail import AuthMail
from datetime import datetime
from pathlib import Path
import json
from extractor.manufactured.datadirs import GENERATED_DIR_INPUTS, INPUT_DIR_GNPS_TASKS

def send():
    INPUT_DIR_GNPS_TASKS.mkdir(parents=True, exist_ok=True)
    
    input_mgf = str((GENERATED_DIR_INPUTS / "All GNPS.mgf").resolve())
    input_quant = str((GENERATED_DIR_INPUTS / "Quantification table.csv").resolve())
    auth = AuthMail.from_txt("../../../Auth GNPS.txt")
    title = f"Manufactured case {datetime.now().isoformat()}"
    (gnps_input_mgf, gnps_input_csv) = _upload_to_gnps(auth, input_mgf, input_quant, title)
    task_ids = invoke_all(auth, gnps_input_mgf, gnps_input_csv, title)

    # ids_by_first = {k:[v for kk,v in e.items()] for k,e in d.items()}
    # all_ids = [v for d in ids_by_first.values() for v in d]
    with open(INPUT_DIR_GNPS_TASKS/"Gnps task ids.json", 'w') as f:
        f.write(json.dumps(task_ids))

def send_1():
    INPUT_DIR_GNPS_TASKS.mkdir(parents=True, exist_ok=True)
    
    input_mgf = str((GENERATED_DIR_INPUTS / "01 GNPS old.mgf").resolve())
    input_quant = str((GENERATED_DIR_INPUTS / "01 Quantification table.csv").resolve())
    auth = AuthMail.from_txt("../../../Auth GNPS.txt")
    title = f"Manufactured case {datetime.now().isoformat()}"
    (gnps_input_mgf, gnps_input_csv) = _upload_to_gnps(auth, input_mgf, input_quant, title)
    task_ids = [invoke(auth, gnps_input_mgf, gnps_input_csv, title, 6, 0.02)]

    with open(INPUT_DIR_GNPS_TASKS/"Gnps task ids.json", 'w') as f:
        f.write(json.dumps(task_ids))

def invoke_all(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description):
    task_ids = []
    for peak in [6, 5, 4]:
        for mass_diff in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]:
            task_id = invoke(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, peak, mass_diff)
            task_ids.append(task_id)
    return task_ids

def invoke(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, peak, mass_diff):
    invokeParameters = _get_networking_parameters()
    invokeParameters["MIN_MATCHED_PEAKS_SEARCH"] = str(peak)
    invokeParameters["MAX_SHIFT_MASS"] = str(mass_diff)
    invokeParameters["desc"] = job_description + \
                '_'+str(peak)+'_'+str(mass_diff)
    invokeParameters["quantification_table"] = "d./" + \
                auth.username + path_file_quan_in_gnps
    invokeParameters["spec_on_server"] = "d./" + \
                auth.username + path_file_mgf_in_gnps
    invokeParameters["email"] = auth.mail
    task_id = _invoke_workflow(auth, "gnps.ucsd.edu", invokeParameters)
    return task_id
