from django.shortcuts import render


def index(request):
    return render(request, 'genbank_fasta/index.html')
