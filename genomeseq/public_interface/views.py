from django.shortcuts import render

# Create your views here.
def index(request):
    return render(request, 'public_interface/index.html')

def browse(request):
    return render(request, 'public_interface/browse.html')
