import itertools

from django.shortcuts import render
from django.db.models import Prefetch

from .models import Vouchers
from .models import FlickrImages


def index(request):
    return render(request, 'public_interface/index.html')


def browse(request):
    queryset = Vouchers.objects.order_by('-timestamp')[:10]

    # TODO improve this ugly hack. Use select_related or prefetch_related
    vouchers_with_images = []
    for i in queryset:
        q = FlickrImages.objects.filter(voucher=i.code)
        if q.count() > 0:
            vouchers_with_images.append(i.code)

    return render(request, 'public_interface/browse.html',
                  {
                      'results': queryset,
                      'vouchers_with_images': vouchers_with_images,
                  },
                  )


def show_voucher(request, voucher_code):
    images_queryset = FlickrImages.objects.filter(voucher=voucher_code)

    queryset = Vouchers.objects.get(code=voucher_code)
    return render(request, 'public_interface/show_voucher.html',
                  {'item': queryset,
                   'images': images_queryset,
                   },
                  )
