import itertools
from django.conf import settings
import datetime
import os
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.views.generic import FormView, DetailView
from django.http import JsonResponse


from .forms import MultiURLForm

class SignUpView(FormView):
    template_name = 'signup.html'
    form_class = MultiURLForm

def test_ajax(request):
    test_string = request.GET.get('test_string', None)
    data = {
        'output': len(test_string)
    }
    return JsonResponse(data)

def index(request):
    return render(request, 'index.html', {'description': "z"})

