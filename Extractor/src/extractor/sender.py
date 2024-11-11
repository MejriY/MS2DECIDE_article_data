from extractor.mgf import Mgf
from ms2decide.AuthMail import AuthMail
from datetime import datetime
from ms2decide.ClosestGNPS import _invoke_workflow
from ms2decide.ClosestGNPS import _get_networking_parameters
from ms2decide.ClosestGNPS import _upload_to_gnps

def send(mgf_file, quantification_file, title_prefix):
    input_mgf = str(mgf_file)
    input_quant = str(quantification_file)
    auth = AuthMail.from_txt("../../../Auth GNPS.txt")
    title = f"{title_prefix} {datetime.now().isoformat()}"
    (gnps_input_mgf, gnps_input_csv) = _upload_to_gnps(auth, input_mgf, input_quant, title)
    return invoke_all(auth, gnps_input_mgf, gnps_input_csv, title)

def invoke_all(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description):
    task_ids = []
    for peak in [6, 5, 4]:
        for mass_diff in [0.02, 0.1, 10, 25, 50, 100, 250, 500, 0]:
            task_id = invoke(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, peak, mass_diff)
            task_ids.append(task_id)
    return task_ids

def invoke(auth, path_file_mgf_in_gnps, path_file_quan_in_gnps, job_description, peak = 6, mass_diff = 0.02):
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
