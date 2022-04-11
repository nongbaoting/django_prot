import subprocess,json
from os import path
import pickle, datetime
from protein.models import SubmitInfoNew
import numpy as np 

def create_tmpDir(baseDir, myuuid):
    tmpdir = path.join(baseDir, myuuid)
    subprocess.run("mkdir -p " + tmpdir, shell=True)
    return tmpdir

def pickle_dump2file(obj, filename): 
    fo =  open(filename, 'wb')
    pickle.dump(list(obj), fo)
    fo.close()

def numpy_save2file(obj, filename):
    np.savez_compressed(filename, cdd=obj)
    


def create_submit_form(params, res):
    # params['uuid'] = res.tasks_id
    # params['task_status'] = res.status
    # params['proj_type'] = res.task_name
    SubmitInfoNew.objects.create(
        uuid = res.task_id,
        proj_type =  params.get('proj_type',default=''),
        email_addr = params.get("email", default=''),
        job_name = params.get("job_name", default='default'),
        upload_date = datetime.datetime.now(),
        tools= params.get('tools', default=''),
        params=params,
        task_status=  res.status
    )

def load_POST(request):
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    return body

