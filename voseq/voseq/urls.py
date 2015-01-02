from django.conf.urls import patterns, include, url
from django.contrib import admin

urlpatterns = patterns(
    '',
    url(r'^blast_local', include('blast_local.urls', namespace='blast_local')),
    url(r'^', include('public_interface.urls', namespace='public_interace')),
    # Examples:
    # url(r'^blog/', include('blog.urls')),

    url(r'^admin/', include(admin.site.urls)),

)
