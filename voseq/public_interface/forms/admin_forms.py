from django.forms import ModelForm, TextInput, ModelChoiceField

from public_interface.models import Vouchers, Sequences, Genes



class VoucherForm(ModelForm):
    code = ModelChoiceField(queryset=Vouchers.objects.distinct("code"),
                            empty_label="Choose a value")

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


class SequenceForm(ModelForm):
    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user", None)
        super(SequenceForm, self).__init__(*args, **kwargs)
        self.fields["gene_code"] = ModelChoiceField(
            queryset=Genes.objects.filter(user=user), empty_label="Choose a value")

    class Meta:
        model =  Sequences
        fields = "__all__"
        exclude = ["user"]


class GeneForm(ModelForm):
    class Meta:
        model = Genes
        fields = "__all__"
        exclude = ["user", "time_created"]
