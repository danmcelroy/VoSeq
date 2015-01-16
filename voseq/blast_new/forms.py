from django import forms


class BLASTNewForm(forms.Form):
    name = forms.CharField(
        label='name',
        max_length=100,
    )
