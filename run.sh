# conda activate webProt
### 更新
# cd /dat1/nbt2/server/PROTsim/public/
# mv django_prot bak/
# cp -r ../django_prot .
cd /dat1/nbt2/server/PROTsim/public/django_prot
conda activate webProt
python run_predict.py watch_blast > logs/blast.log 2>&1 &
celery -A django_prot worker -l info  -c 4 -E > logs/celery.log 2>&1 &
python manage.py runserver 222.200.186.124:3389


# dev
conda activate webProt
cd /dat1/nbt2/server/PROTsim/django_prot
python manage.py runserver 222.200.186.124:9002