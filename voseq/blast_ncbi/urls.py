from django.conf.urls import patterns
from django.conf.urls import url

from . import views


urlpatterns = patterns(
    '',
    url(r'^(?P<voucher_code>.+)/(?P<gene_code>.+)/$', views.index, name='index'),
)
