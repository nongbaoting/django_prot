from django.urls import path
from .views import browse, superposed, struc_align
from .views import browse_cas12f1
app_name = "crispr"
urlpatterns = [
    path("api/browse", browse.browse),
    path("api/structure/getFile/", browse.struc_getFile),#cas9f
    path("api/structure/alignScore/TMScore/", browse.alignTMscore),
    path("api/structure/alignScore/SPScore/", browse.alignSPscore),
    path("api/structure/alignScore/FatcatScore/", browse.alignFatcatScore),
    path("api/structure/cas12f1/FatcatScore/", browse_cas12f1.alignFatcatScore), #cas12f
    path("api/structure/superposed/pairs/",   superposed.alignment),
    path("api/similarity/aligment/TMalign", struc_align.pair_TMalign),
    path("api/similarity/aligment/SPalign", struc_align.pair_SPalign),
    path("api/similarity/aligment/Fatcat", struc_align.pair_Fatcat),
]
