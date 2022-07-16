from django.urls import path
from .views import browse

app_name = "crispr"
urlpatterns = [
    path("api/browse", browse.browse),
    path("api/structure/getFile/", browse.struc_getFile),
    path("api/structure/alignScore/TMScore/", browse.alignTMscore),
    path("api/structure/alignScore/SPScore/", browse.alignSPscore),
]
