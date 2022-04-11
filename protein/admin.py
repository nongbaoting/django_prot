from django.contrib import admin
from .models import *
# Register your models here.
# 注册以便 admin网页查看
admin.site.register(SubmitInfo)
admin.site.register(GeneInfo)
admin.site.register(Structure)
admin.site.register(AF_uniprot)
admin.site.register(SubmitInfoNew)
admin.site.register(ProtCDD)