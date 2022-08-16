from django.urls import path
from .views import browse, superposed, struc_align

app_name = "crispr"
urlpatterns = [
    path("api/browse", browse.browse),
    path("api/structure/getFile/", browse.struc_getFile),
    path("api/structure/alignScore/TMScore/", browse.alignTMscore),
    path("api/structure/alignScore/SPScore/", browse.alignSPscore),
    path("api/structure/alignScore/FatcatScore/", browse.alignFatcatScore),
    path("api/structure/superposed/pairs/",   superposed.alignment),
    path("api/similarity/aligment/TMalign", struc_align.pair_TMalign),
    path("api/similarity/aligment/SPalign", struc_align.pair_SPalign),
    path("api/similarity/aligment/Fatcat", struc_align.pair_Fatcat),
]
