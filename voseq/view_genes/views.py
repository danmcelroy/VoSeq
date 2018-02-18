from django.shortcuts import render

from core.utils import get_context
from public_interface.models import Genes
from stats.models import VouchersPerGene


def index(request):
    context = get_context(request)
    voucher_count = get_voucher_count()

    queryset = Genes.objects.all().values()
    result = []
    for i in queryset:
        gene_code = i['gene_code']
        try:
            i['voucher_count'] = voucher_count[gene_code]
        except KeyError:
            continue
        result.append(i)
    result = sorted(result, key=lambda k: k['gene_code'].lower())

    context["result"] = result
    return render(request, 'view_genes/index.html', context)


def get_voucher_count():
    voucher_count = {}
    queryset = VouchersPerGene.objects.all().values('gene_code', 'voucher_count')
    for i in queryset:
        gene_code = i['gene_code']
        voucher_count[gene_code] = i['voucher_count']
    return voucher_count


def gene(request, gene_code):
    context = get_context(request)

    queryset = Genes.objects.filter(gene_code=gene_code)
    if not queryset:
        item = ''
    else:
        item = queryset[0]

        if ';' in item.intron:
            introns = item.intron.split(';')
            j = 1
            out = ''
            for i in introns:
                out += '<b>Intron ' + str(j) + ':</b> ' + i + ' > '
                j += 1
            item.intron = out.rstrip(' > ')

    context["item"] = item
    return render(request, 'view_genes/gene.html', context)
