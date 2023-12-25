from django.urls import path
from .views import struc_similarity
from .views import struc_predict, struc_align
from .views import blast
from .views import login
from .views import cdd, cdd_search
from .views import phylogenetic
from .views import jobQueue
from .views import domain_pdbe
from .views import pdb_domain_annotations
from .views.results import tada_like,csr, repeatDomain
from .views import test

app_name = 'protein'
urlpatterns = [
    path('login', login.login),

    # Structure predict
    path('get_token', struc_predict.get_token, name="get_token"),
    path('', struc_predict.index, name='index'),
    path('predict/structure_submit/', struc_predict.structure_submit, name="structure_submit"),
     
    path('predict/structure/result',
         struc_predict.structure_result, name="structure_result"),
    path('predict/structure/result', struc_predict.alphafold2, name="alphafold2"),
    path('predict/structure/check/', struc_predict.check_proj, name="check_proj"),
    path('api/structure/queue/', struc_predict.struc_queue, name='struc_queue'),
    path('api/structure/search/', struc_predict.struc_search, name='struc_search'),
    path("api/structure/getFile/",
         struc_predict.struc_getFile, name='struct_getFile'),
    path("api/structure/getTemplate/", struc_predict.struc_template),

    # new
     path('predict/structure_submit_new/', struc_predict.submit_new),
     path("api/structure/getFile_new/",    struc_predict.getFile),
     path('api/structure/result_new/',  struc_predict.get_result),
      path("api/structure/getTemplate_new/", struc_predict.struc_template_new),

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
    path("api/similarity/aligment/TMalign", struc_align.pair_TMalign),
    path("api/similarity/aligment/SPalign", struc_align.pair_SPalign),
    path("api/similarity/aligment/Fatcat", struc_align.pair_Fatcat),

    # sequnce similarity
    path("api/blast/psijackhmmer/", blast.psijackhmmer),
    path("api/blast/res/blast_jackhmmer/", blast.res_blast_jackhmmer),
     path("api/blast/res/architecture/", blast.architecture),
    path("api/blast/queue/", blast.queue),
    path("api/blast/queue/search/", blast.search),

    # cdd
    path("api/cdd/search/", cdd.search),
    path("api/cdd/search/page/", cdd.pages),
    path("api/cdd/search/get_all_protin_ids/", cdd.get_all_protin_id),
    path("api/cdd/search_save/", cdd_search.search_save),
    path("api/cdd/retrieve_save/", cdd_search.retrieve_save),
    path("api/cdd/filter_cdd_save/", cdd_search.filter_cdd_save),

    # phylogenectics
    path("api/phylo/run/", phylogenetic.run_phylo),
    path("api/phylo/get_fasta/", phylogenetic.get_fasta),
    path("api/phylo/tree_files/", phylogenetic.tree_files),
    # queue
    path('api/Queue/', jobQueue.jobTable),

    # domain_annotation
    path("api/pdbe/interpro/",domain_pdbe.interpro),
    path("api/pdb_domain_annotations/uploadPDB_and_annotation/", pdb_domain_annotations.uploadPDB_and_annotation),
    path("api/pdb_domain_annotations/parser_results/", pdb_domain_annotations.parser_results),
    path("api/pdb_domain_annotations/get_pdbFile/", pdb_domain_annotations.get_pdbFile),
    path("api/pdb_domain_annotations/align/", pdb_domain_annotations.align),

    # tadA_like, results
    path("api/results/tada_like/", tada_like.interpro),
    path("api/results/get_pdbFile/", tada_like.get_pdbFile),
     path("api/results/csr_like/", csr.interpro),
    path("api/results/csr_pdbFile/", csr.get_pdbFile),
    ## repeat proteins
     path("api/results/repeatDomain/fetchData/", repeatDomain.fetchData),
     path("api/results/repeatDomain/fetchProteins/", repeatDomain.fetchProteins),


    # test
    path("api/test/get_pdbFile/", test.get_pdbFile),


   

]
