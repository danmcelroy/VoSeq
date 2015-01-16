from django import forms


class BLASTNewForm(forms.Form):
    name = forms.CharField(
        label='name',
        max_length=100,
    )
    sequence = forms.CharField(
        label='sequence',
        required=True,
        widget=forms.Textarea,
    )
