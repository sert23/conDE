from django import forms


class MultiURLForm(forms.Form):
    SRRtext = forms.CharField(label="Paste SRA IDs (starting with SRR or ERR, one per line) ", widget=forms.Textarea, required=False)
    URLtext = forms.CharField(label="Paste URL/links with the files (one per line) ", widget=forms.Textarea, required=False)

