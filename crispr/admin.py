from django.contrib import admin

# Register your models here.
from .models import *

admin.site.register(CASInfo)
admin.site.register(AlignSPScore)
admin.site.register(AlignTMScore)
