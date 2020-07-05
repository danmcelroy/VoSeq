from django.urls import path

from api import views

urlpatterns = [
    path('simple-search/', views.simple_search)
]
