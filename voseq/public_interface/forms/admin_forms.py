from django.forms import ModelForm, TextInput

from public_interface.models import Vouchers



class VoucherForm(ModelForm):
    class Meta:
        model = Vouchers
        fields = "__all__"
        exclude = ["edits", "latest_editor", "user"]
        widgets = {
            "orden": TextInput,
            "superfamily": TextInput,
            "family": TextInput,
            "subfamily": TextInput,
            "tribe": TextInput,
            "subtribe": TextInput,
            "genus": TextInput,
            "species": TextInput,
            "subspecies": TextInput,
            "country": TextInput,
            "specific_locality": TextInput,
            "collector": TextInput,
            "extraction": TextInput,
            "extraction_tube": TextInput,
            "extractor": TextInput,
            "voucher_locality": TextInput,
            "published_in": TextInput,
            "notes": TextInput,
            "hostorg": TextInput,
            "voucher_code": TextInput,
            "code_bold": TextInput,
            "determined_by": TextInput,
            "author": TextInput,
        }
