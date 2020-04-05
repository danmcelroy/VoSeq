from django.conf.urls import url

from . import views

app_name = 'create_dataset'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^results/$', views.results, name='results'),
    url(r'^results/(?P<file_name>.+\.txt)/$', views.serve_file, name='serve_file'),
]
