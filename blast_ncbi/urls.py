from django.conf.urls import url

from . import views


app_name = 'blast_local_ncbi'
urlpatterns = [
    url(r'^(?P<voucher_code>.+)/(?P<gene_code>.+)/$', views.index, name='index'),
]
