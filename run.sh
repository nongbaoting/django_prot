conda activate webProt
celery -A django_prot worker -l info  -c 4 -E > logs/celery.log 2>&1 &
python manage.py runserver 222.200.186.124:3389