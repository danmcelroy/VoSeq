from django.forms import ModelForm

from public_interface.models import Vouchers



class VoucherForm(ModelForm):
    class Meta:
        model = Vouchers
        fields = ["code"]
