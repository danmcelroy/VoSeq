from django.conf.urls import patterns
from django.conf.urls import url
from django.conf.urls import include

from haystack.query import SearchQuerySet
from haystack.views import SearchView
from haystack.views import search_view_factory

from . import views
from .forms import AdvancedSearchForm


urlpatterns = patterns(
    '',
    url(r'^$', views.index, name='index'),
    url(r'^browse/$', views.browse, name='browse'),
    url(r'^search/$', views.search, name='search'),
    # url(r'^search/', include('haystack.urls')),
    url(r'^p/(?P<voucher_code>.+)/$', views.show_voucher, name='show_voucher'),
    url(r'^s/(?P<voucher_code>.+)/(?P<gene_code>.+)/$', views.show_sequence, name='show_sequence'),
)

urlpatterns += patterns(
    'haystack.views',
    url(r'^$', SearchView(
        template='search/search.html',
        form_class=AdvancedSearchForm,
    ), name='haystack_search'),
)
