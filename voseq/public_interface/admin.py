from django.contrib import admin
from django.core.exceptions import PermissionDenied
from django.http import HttpRequest

from public_interface.models import Vouchers
from public_interface.views import change_selected


# Customize what and the way you show it
class VouchersAdmin(admin.ModelAdmin):
    list_display = ['code', 'genus', 'species', 'sex', 'voucher', 'country', 'collector']
    ordering = ['code']
    list_filter = ['voucher']
    search_fields = ['=genus', '=species']

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

                 (None, {'fields': ['timestamp']}),
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


# Register your models here.
admin.site.register(Vouchers, VouchersAdmin)
