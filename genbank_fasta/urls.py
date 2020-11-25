from django.conf.urls import url

from . import views


app_name = 'genbank_fasta'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^results/$', views.generate_results, name='generate-genbank-fasta-results'),
    url(r'^results/(?P<dataset_id>[0-9]+)/$', views.results, name='create-genbank-results'),
    url(r'^download/(?P<dataset_id>[0-9]+)/$', views.serve_file, name='download-genbank-fasta-results'),
]
