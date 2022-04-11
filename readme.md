```bash
# django-admin startproject mysite

# 创建应用
python manage.py startapp protein

# 应用数据库
# del 
from protein.models import *; all = SubmitInfoNew.objects.all() ;  all.delete()
all = SubmitInfo.objects.all()
all.delete()

# 更改
python manage.py makemigrations protein
python manage.py migrate

python txt2db.py geneinfo  ../django_21/sprot_info.txt
python txt2db.py Structure ../django_21/sprot_file_relate.txt

# 启用服务
pip install livereload -i http://mirrors.aliyun.com/pypi/simple/ --trusted-host mirrors.aliyun.com

python manage.py runserver 222.200.186.47:9000

# publish test
python manage.py runserver 222.200.186.47:9001

# 启动

# 1. 启动 java. Biozernike
cp -r /training/nong/web/Dev/django_prot .
cd /training/nong/web/java
nohub java -cp /training/nong/web/java/CalProSimilariry-1.0.1-jar-with-dependencies.jar sysu.JPype.Compare

# 2. 启动wachdog。 alphafold, roseTTAFold
cd /training/nong/web/public/django_prot
source activate webProt
python ./run_predict.py watchdog > logs/logs.txt 2>&1 &
# 2.1 启动blast 124
cd /training/nong/web/public/django_prot
source activate webProt
python run_predict.py  watch_blast /dat1/nbt2/proj/21-prot/web/data/uploads/blast

# 3. 启动 Django
sshfs -o nonempty nbt2@222.200.186.124:/dat1/ /dat1/
source activate webProt
cd /training/nong/web/public/django_prot
python manage.py runserver 222.200.186.47:9000

# apache
sudo systemctl restart apache2

# dev
python manage.py runserver 222.200.186.47:9001
```

# django 目录
`/training/nong/web/public`

# dev to public
cd /training/nong/web/public/django_prot
cp db.sqlite3 .. 
cp -r  /training/nong/web/Dev/django_prot/* .
cp ../db.sqlite3 .
python manage.py makemigrations protein
python manage.py migrate

rm -rf /training/nong/protein/db/web/protein
cp -r /training/nong/web/Dev/prot/dist/ /training/nong/protein/db/web/protein

# 2.1 启动blast 124
cd /training/nong/web/public/django_prot
source activate webProt
python run_predict.py  watch_blast 

#/dat1/nbt2/proj/21-prot/web/data/uploads/blast


# mysql
mysqldump -u nongbt -pNBT9175.814@lys628 --no-tablespaces --databases nongbt_db > /dat1/nbt2/proj/21-prot/web/data/mysql/nongbt_db`dateymd`.sql

# mysql
drop database nongbt_db;
CREATE DATABASE nongbt_db CHARACTER SET utf8 COLLATE utf8_general_ci;
# load
mysql -u nongbt -pNBT9175.814@lys628  < /dat1/nbt2/proj/21-prot/web/data/mysql/nongbt_db220326.sql
# source /home/abc/abc.sql 
# rm protein/migrations/*