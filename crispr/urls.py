from django.urls import path
from .views import browse, superposed

app_name = "crispr"
urlpatterns = [
    path("api/browse", browse.browse),
    path("api/structure/getFile/", browse.struc_getFile),
    path("api/structure/alignScore/TMScore/", browse.alignTMscore),
    path("api/structure/alignScore/SPScore/", browse.alignSPscore),
    path("api/structure/alignScore/FatcatScore/", browse.alignFatcatScore),
    path("api/structure/superposed/pairs/",   superposed.alignment),
]
