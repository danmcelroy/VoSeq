from django.conf.urls import url

from . import views


urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^(?P<gene_code>.+)/$', views.gene, name='gene'),
]
