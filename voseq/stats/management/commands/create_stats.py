from django.core.management.base import BaseCommand
from django.core.management.base import CommandError

from public_interface.models import Voucher
from public_interface.models import Sequence
from stats.models import Stats


class Command(BaseCommand):
    help = 'Extracts total number of orders, families, genera, etc. from ' \
           'our database.'

    def handle(self, *args, **options):
        queryset = Voucher.objects.all()

        orders = set()
        families = set()
        genera = set()
        species = set()

        for i in queryset:
            order = i.orden
            family = i.family
            genus = i.genus
            this_species = i.genus + ' ' + i.species

            if order != '':
                orders.add(order)
            if family != '':
                families.add(family)
            if genus != '':
                genera.add(genus)
            if this_species != '':
                species.add(this_species)

        num_vouchers = len(queryset)
        num_orders = len(orders)
        num_families = len(families)
        num_genera = len(genera)
        num_species = len(species)

        queryset = Sequence.objects.all()
        num_sequences = len(queryset)

        Stats.objects.update_or_create(
            id=1,
            defaults={
                     'vouchers': num_vouchers,
                     'orders': num_orders,
                     'families': num_families,
                     'genera': num_genera,
                     'species': num_species,
                     'sequences': num_sequences,
                     }
        )
