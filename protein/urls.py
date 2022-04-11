from operator import imod
from django.urls import path
from .views import struc_similarity
from .views import struc_predict
from .views import blast
from .views import login
from .views import cdd, cdd_search
from .views import phylogenetic
from .views import jobQueue

app_name = 'protein'
urlpatterns = [
     path('login',login.login ),
    # Structure predict
    path('get_token', struc_predict.get_token, name="get_token"),
    path('', struc_predict.index, name='index'),
    path('predict/structure_submit/',
         struc_predict.structure_submit, name="structure_submit"),
    path('predict/structure/result',
         struc_predict.structure_result, name="structure_result"),
    path('predict/structure/result', struc_predict.alphafold2, name="alphafold2"),
    path('predict/structure/check/', struc_predict.check_proj, name="check_proj"),
    path('api/structure/queue/', struc_predict.struc_queue, name='struc_queue'),
    path('api/structure/search/', struc_predict.struc_search, name='struc_search'),
    # search
    path('search/', struc_predict.search, name="search"),
    path("test/", struc_predict.test, name="test"),
    path("show_structure/", struc_predict.show_structure, name="show_structure"),
    #path("show_structure/<str:str_id>/", views.show_structure, name="show_str_detail")

    # Structure similarity
    path("api/similarity/upload_pdb/", struc_similarity.upload_pdb), 
    path("api/similarity/DUF_SPalign/", struc_similarity.DUF_SPalign), 
    path("api/similarity/results/", struc_similarity.results),
    path("api/similarity/getOneItem/", struc_similarity.getOneItem),

    # sequnce similarity
    path("api/blast/psijackhmmer/", blast.psijackhmmer),
    path("api/blast/res/blast_jackhmmer/", blast.res_blast_jackhmmer),
    path("api/blast/queue/", blast.queue),
    path("api/blast/queue/search/", blast.search),
    
    # cdd 
    path("api/cdd/search/", cdd.search),
    path("api/cdd/search/page/", cdd.pages),
    path("api/cdd/search/get_all_protin_ids/", cdd.get_all_protin_id),
    path("api/cdd/search_save/", cdd_search.search_save), 
    path("api/cdd/retrieve_save/", cdd_search.retrieve_save), 
    path("api/cdd/filter_cdd_save/", cdd_search.filter_cdd_save), 
    #phylogenectics
     path("api/phylo/run/", phylogenetic.run_phylo),
     path("api/phylo/get_fasta/", phylogenetic.get_fasta),
     path("api/phylo/tree_files/", phylogenetic.tree_files), 
     #queue
     path('api/Queue/', jobQueue.jobTable),
]
