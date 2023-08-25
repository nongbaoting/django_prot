# 创建

django-admin startproject mysite

# 创建应用

python manage.py startapp protein

# 应用数据库

## 更改

python manage.py makemigrations protein
python manage.py migrate

python txt2db.py geneinfo  ../django_21/sprot_info.txt
python txt2db.py Structure ../django_21/sprot_file_relate.txt

## del

from protein.models import *; all = SubmitInfoNew.objects.all() ;  all.delete()
all = SubmitInfo.objects.all()
all.delete()

# 启用服务

pip install livereload -i http://mirrors.aliyun.com/pypi/simple/ --trusted-host mirrors.aliyun.com

python manage.py runserver 222.200.186.47:9000

# publish test

python manage.py runserver 222.200.186.47:9001

#### 重现步骤 安装

conda env create -f environment.yml

# 启动

发布 ，更改redis `django_prot/celeryconfig.py`

# 1. 启动 java. Biozernike

<!-- cp -r /training/nong/web/Dev/django_prot .
cd /training/nong/web/java
nohub java -cp /training/nong/web/java/CalProSimilariry-1.0.1-jar-with-dependencies.jar sysu.JPype.Compare -->

# 2. 启动wachdog。 alphafold, roseTTAFold

cd /training/nong/web/public/django_prot
source activate webProt
python ./run_predict.py watchdog > logs/logs.txt 2>&1 &

# 2.1 启动blast 124

cd /training/nong/web/public/django_prot
source activate webProt
python run_predict.py  watch_blast /dat1/nbt2/proj/21-prot/web/data/uploads/blast

#### 重现步骤 3. 启动 Django

sshfs -o nonempty nbt2@222.200.186.124:/dat1/ /dat1/
source activate webProt
cd /training/nong/web/public/django_prot
python manage.py runserver 222.200.186.47:9000

# apache

sudo systemctl restart apache2

# dev

python manage.py runserver 222.200.186.47:9001

# dev to public

cd /dat1/nbt2/server/PROTsim/public

mv django_prot bak/
cp -r ../django_prot .

# 2.1 启动blast 124

cd /training/nong/web/public/django_prot
source activate webProt
python run_predict.py  watch_blast

#/dat1/nbt2/proj/21-prot/web/data/uploads/blast

# mysql

mysqldump -u nongbt -pNBT9175.814@lys628 --no-tablespaces --databases nongbt_db > /dat1/nbt2/proj/21-prot/web/data/mysql/nongbt_db `dateymd`.sql

drop database nongbt_db;
CREATE DATABASE nongbt_db CHARACTER SET utf8 COLLATE utf8_general_ci;

# load

mysql -u nongbt -pNBT9175.814@lys628  < /dat1/nbt2/proj/21-prot/web/data/mysql/nongbt_db220326.sql

# source /home/abc/abc.sql

# rm protein/migrations/*


# structure annotation

## binding site
* GeoBind ,structure-base (server 150) 

    https://github.com/zpliulab/GeoBind

* clape # protein-DNA RNA ligand , #sequence-base

    https://github.com/YAndrewL/CLAPE/tree/main
    python clape.py --input example.fa --output out.txt

* scannet #* PPI, #structure-base (server 150)

* mymetal metal ！！！！#sequence-base

    https://academic.oup.com/bioinformatics/article/38/14/3532/6594112?login=false#401869179

* LMetalSite metals #sequence-base
    https://github.com/biomed-AI/LMetalSite