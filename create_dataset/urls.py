from django.conf.urls import url

from . import views

app_name = 'create_dataset'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^results/$', views.generate_results, name='generate-dataset-results'),
    url(r'^results/(?P<dataset_id>[0-9]+)/$', views.results, name='create-dataset-results'),
    url(r'^download/(?P<dataset_id>[0-9]+)/$', views.serve_file, name='download-dataset-results'),
]
