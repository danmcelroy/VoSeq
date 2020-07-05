from django.conf.urls import url

from . import views


app_name = 'overview_table'
urlpatterns = [
    url(r'^$', views.index, name='index'),
]
