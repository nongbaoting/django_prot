import subprocess,json
from os import path
import pickle, datetime
from protein.models import SubmitInfoNew
import numpy as np 
from collections import defaultdict

def create_tmpDir(baseDir, myuuid):
    tmpdir = path.join(baseDir, myuuid)
    subprocess.run("mkdir -p " + tmpdir, shell=True)
    return tmpdir

def pickle_dump2file(obj, filename): 
    fo =  open(filename, 'wb')
    pickle.dump(list(obj), fo)
    fo.close()

def pickle_dumpObj2file(obj, filename): 
    fo =  open(filename, 'wb')
    pickle.dump(obj, fo)
    fo.close()

def pickle_load_file(filename):
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj

def numpy_save2file(obj, filename):
    np.savez_compressed(filename, cdd =obj)

def run_cmd(cmd):
    std = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return std

def create_submit_form(params, res):
    # params['uuid'] = res.tasks_id
    # params['task_status'] = res.status
    # params['proj_type'] = res.task_name
    SubmitInfoNew.objects.create(
        uuid = res.task_id,
        proj_type =  params.get('proj_type',''),
        email_addr = params.get("email", ''),
        job_name = params.get("job_name",'default'),
        upload_date = datetime.datetime.now(),
        tools= params.get('tools', ''),
        params=params,
        task_status=  res.status
    )

def load_POST(request):
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    return body

colorSet = [
  '#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462', '#B3DE69',
  '#FCCDE5',
  '#D9D9D9',
   '#D9D9D9',
  '#CCEBC5',
  '#FFED6F',
  '#1B9E77',
  '#D95F02',
  '#7570B3',
  '#E7298A',
  '#66A61E',
  '#E6AB02',
  '#A6761D',
  '#666666',
  '#543005',
  '#8C510A',
  '#BF812D',
  '#DFC27D',
  '#F6E8C3',
  '#C7EAE5',
  '#80CDC1',
  '#35978F',
  '#01665E',
  '#003C30',
]

colorSet3D = [
    "#27A3B4",
    # "#C08423",
    "#D83E7C",
    "#1B9E77",


'#FCCDE5',
  '#D9D9D9',
  '#BC80BD',
  '#CCEBC5',
  '#FFED6F',
  '#1B9E77',
  '#D95F02',
  '#7570B3',
  '#E7298A',
  '#66A61E',
  '#E6AB02',
  '#A6761D',
]

def addCDcorlor(content,colors=colorSet):
    myd = defaultdict(str)
    index = 0
    new_arr = []
    for items in content:
        # print(items['cdd_locs'])
        for loc in items['cdd_locs']:
            if not loc[2] in myd:
                myd[loc[2]] =  colors[index]
                index += 1
                if index >= len(colors):
                    index=0
            loc.append(myd[ loc[2] ])
        new_arr.append(items)
    return new_arr

def addCDArchiColor(content,colors=colorSet):
    myd = defaultdict(str)
    index = 0
    new_arr = []
    for item in content:
        cddNames = item["cdd_nameCat"].split(', ')
        cd_color = []
        for cd_name in cddNames:
            if cd_name not in myd:
                myd[cd_name] = colors[index]
                index += 1
                if index >= len(colors):
                        index=0
            cd_color.append(myd[cd_name])
        item['colors'] = cd_color
        item['cddNames'] = cddNames
        new_arr.append(item)
    
    return new_arr
