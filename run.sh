# conda activate webProt
# /dat1/nbt2/server/PROTsim/public/django_prot
python run_predict.py watch_blast > logs/blast.log 2>&1 &
celery -A django_prot worker -l info  -c 4 -E > logs/celery.log 2>&1 &
python manage.py runserver 222.200.186.124:3389