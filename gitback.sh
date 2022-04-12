#/bin/bash
cd /dat1/nbt2/server/PROTsim/django_prot
git add -A
git commit -m 'routine backup '`date +%Y%m%d`
git push -u origin main
cd /dat1/nbt2/server/PROTsim/prot
git add -A
git commit -m 'routine backup '`date +%Y%m%d`
git push -u origin main
