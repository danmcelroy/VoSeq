from django.conf.urls import patterns, include, url
from django.contrib import admin

urlpatterns = patterns(
    '',
    url(r'^blast_local', include('blast_local.urls', namespace='blast_local')),
    url(r'^blast_local_full', include('blast_local_full.urls', namespace='blast_local_full')),
    url(r'^', include('public_interface.urls', namespace='public_interace')),

    url(r'^admin/', include(admin.site.urls)),

)
