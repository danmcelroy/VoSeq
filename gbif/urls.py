from django.conf.urls import url

from . import views


app_name = 'gbif'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^dump_data/$', views.dump_data, name='dump_data'),
]
