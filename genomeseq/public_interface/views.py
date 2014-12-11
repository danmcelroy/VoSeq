from django.shortcuts import render

from .models import Vouchers
from .models import FlickrImages


def index(request):
    return render(request, 'public_interface/index.html')


def browse(request):
    queryset = Vouchers.objects.all().order_by('-timestamp')[:10]
    return render(request, 'public_interface/browse.html',
                  {'results': queryset},
                  )


def show_voucher(request, voucher_code):
    images_queryset = FlickrImages.objects.filter(voucher=voucher_code)

    queryset = Vouchers.objects.get(code=voucher_code)
    return render(request, 'public_interface/show_voucher.html',
                  {'item': queryset,
                   'images': images_queryset,
                   },
                  )
