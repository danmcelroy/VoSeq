from django import forms
from django.conf import settings
from django.contrib import admin
from django.core.exceptions import PermissionDenied
from django.db import models
from django.http import HttpRequest

from public_interface.models import FlickrImages
from public_interface.models import LocalImages
from public_interface.models import TaxonSets
from public_interface.models import Sequences
from public_interface.models import Vouchers
from public_interface.views import change_selected


class ImageInLine(admin.StackedInline):
    model = LocalImages
    fields = ['voucherImage']


class FlickImageInLine(admin.StackedInline):
    model = FlickrImages
    fields = ['voucherImage']


# Customize what and the way you show it
class VouchersAdmin(admin.ModelAdmin):
    list_display = ['code', 'genus', 'species', 'sex', 'voucher', 'country', 'collector']
    ordering = ['code']
    search_fields = ['=code', '=genus', '=species']

    # list_editable = ['genus', 'species', 'sex', 'voucher', 'country', 'collector']

    actions = ['batch_changes']

    fieldsets = [('Voucher Information', {'fields': ['code', 'voucher', 'voucherLocality',
                                                     'voucherCode']}
                  ),

                 ('Specimen Information', {'fields': ['orden', 'superfamily', 'family',
                                                      'subfamily', 'tribe', 'subtribe',
                                                      'genus', 'species', 'subspecies',
                                                      'hostorg', 'author', 'typeSpecies',
                                                      ],
                                           'classes': ['collapse']}),

                 ('Collection Information', {'fields': ['country', 'specificLocality',
                                                        'latitude', 'longitude',
                                                        'max_altitude', 'min_altitude',
                                                        'collector', 'code_bold',
                                                        'dateCollection', 'determinedBy',
                                                        'sex', 'extractor', 'extraction',
                                                        'extractionTube', 'notes',
                                                        'publishedIn', 'dateExtraction',
                                                        'edits', 'latesteditor',
                                                        ],
                                             'classes': ['collapse']}),
                 ]

    def batch_changes(self, request, queryset):
        # Check that the user has change permission for the actual model
        if not self.has_change_permission(request):
            raise PermissionDenied
        else:
            selected = request.POST.getlist(admin.ACTION_CHECKBOX_NAME)
            new_request = HttpRequest()
            new_request.method = 'GET'
            return change_selected(new_request, ",".join(selected))

    batch_changes.short_description = "Change selected in batch"
    formfield_overrides = {
        models.TextField: {'widget': forms.TextInput}
    }

    inlines = []

    if settings.PHOTOS_REPOSITORY == 'flickr':
        inlines.append(FlickImageInLine)
    else:
        inlines.append(ImageInLine)


class SequencesAdmin(admin.ModelAdmin):
    # TODO let users know that code and genecode keywords act as AND boolean search
    search_fields = ['=code__code', '=gene_code']
    list_display = ['code', 'gene_code', 'genbank', 'accession', 'labPerson', 'notes', 'time_edited', 'time_created']
    fields = ['code', 'gene_code', 'sequences', 'genbank', 'accession', 'labPerson', 'notes']


class TaxonSetsAdmin(admin.ModelAdmin):
    list_display = ['taxonset_name', 'taxonset_creator', 'taxonset_description']


# Register your models here.
admin.site.register(Sequences, SequencesAdmin)
admin.site.register(TaxonSets, TaxonSetsAdmin)
admin.site.register(Vouchers, VouchersAdmin)
