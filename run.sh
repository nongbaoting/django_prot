#  发布更改
/dat1/nbt2/server/PROTsim/django_prot/django_prot/celeryconfig.py
broker_url ='redis://localhost:6379/1'

# conda activate webProt
### 更新
# cd /dat1/nbt2/server/PROTsim/public/
# mv django_prot bak/
# cp -r ../django_prot .
cd /dat1/nbt2/server/PROTsim/public/django_prot
conda activate webProt
#python run_predict.py watch_blast > logs/blast.log 2>&1 &
celery -A django_prot worker -l info  -c 4 -E > logs/celery.log 2>&1 &
python manage.py runserver 172.22.148.191:3389

#npm
cd /dat1/nbt2/server/PROTsim/prot
npm run build
rm -rf /training/nong/protein/db/web/protein
cp -r /dat1/nbt2/server/PROTsim/prot/dist/ /training/nong/protein/db/web/protein

# 150(47) 启动 alphafold
## cp -r /dat1/nbt2/server/PROTsim/public/django_prot /training/nong/web/django_prot
cd /training/nong/web/django_prot
conda activate webProt2
python ./run_predict.py watchdog > logs/logs.txt 2>&1 &

# dev
conda activate webProt
cd /dat1/nbt2/server/PROTsim/django_prot
python manage.py runserver 172.22.148.191:9002


# 备份数据库
# mysqldump -uroot -p  nongbt_db >nongbt_db.sql
# mysql -uroot -p  nongbt_db <nongbt_db.sql