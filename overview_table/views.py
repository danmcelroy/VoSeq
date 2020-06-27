from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.shortcuts import render

from core.utils import get_context
from public_interface.models import Genes
from overview_table.models import OverviewTable


def index(request):
    context = get_context(request)

    genes = Genes.objects.all().order_by('gene_code')
    vouchers = OverviewTable.objects.all()

    paginator = Paginator(vouchers, 100)

    page = request.GET.get('page')
    try:
        vouchers_for_page = paginator.page(page)
    except PageNotAnInteger:
        vouchers_for_page = paginator.page(1)
    except EmptyPage:
        vouchers_for_page = paginator.page(paginator.num_pages)

    context["data"] = vouchers_for_page
    context["genes"] = genes
    context["page_range"] = paginator.page_range
    return render(request, 'overview_table/index.html', context)
