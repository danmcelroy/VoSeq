from django.shortcuts import render
from django.core.context_processors import csrf
from django.contrib.auth.views import login


def login(request):
    return render(request, 'accounts/login.html')
    # url(r'^accounts/login/$', 'django.contrib.auth.views.login', {'template_name': 'accounts/login.html'}),
