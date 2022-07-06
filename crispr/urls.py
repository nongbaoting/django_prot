from django.urls import path
from .views import browse
app_name = "crispr"
urlpatterns = [
    path("api/browse", browse.browse),
]
