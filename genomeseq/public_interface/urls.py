from django.conf.urls import patterns, url

from . import views

urlpatterns = patterns(
    '',
    url(r'^$', views.index, name='index'),
    url(r'^browse$', views.browse, name='browse'),
    url(r'^p/(?P<voucher_code>.+)/$', views.show_voucher, name='show_voucher'),
)
