import itertools

from django.http import Http404
from django.shortcuts import render

from .models import Vouchers
from .models import FlickrImages
from .models import Sequences


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
    try:
        queryset = Vouchers.objects.get(code__iexact=voucher_code)
    except Vouchers.DoesNotExist:
        raise Http404

    images_queryset = FlickrImages.objects.filter(voucher=voucher_code)

    seqs_queryset = Sequences.objects.filter(code=voucher_code).order_by('gene_code')
    for item in seqs_queryset:
        seq = item.sequences
        item.sequence_length = len(seq)
        item.ambiguous_seq_length = seq.count('?') + seq.count('-') + seq.count('N') + seq.count('n')
        if item.labPerson is not None:
            item.labPerson = item.labPerson.split(" ")[0]

    return render(request, 'public_interface/show_voucher.html',
                  {'item': queryset,  # TODO change item to voucher
                   'images': images_queryset,
                   'sequences': seqs_queryset,
                   },
                  )


def show_sequence(request, voucher_code, gene_code):
    try:
        queryset = Vouchers.objects.get(code__iexact=voucher_code)
    except Vouchers.DoesNotExist:
        raise Http404

    seqs_queryset = Sequences.objects.get(code=voucher_code, gene_code=gene_code)
    images_queryset = FlickrImages.objects.filter(voucher=voucher_code)

    return render(request, 'public_interface/show_sequence.html',
                  {
                      'voucher': queryset,
                      'sequence': seqs_queryset,
                      'images': images_queryset,
                  },)
