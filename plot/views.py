from django.shortcuts import render, redirect
from conDE.settings import MEDIA_ROOT, MEDIA_URL, RSCRIPT_PATH, LOCAL_TEST, RPLOTS_PATH
import pandas as pd
import os
from django.http import JsonResponse
# Create your views here.
from django.views.generic import FormView, DetailView
import json



def random_redirect(request):
    id = request.GET.get('id', None)
    plot = request.GET.get('plot', None)
    redirect(MEDIA_URL,id,plot)

def heatmap_recalculate(request):
    id = request.GET.get('id', None)
    folder = os.path.join(MEDIA_ROOT, id)
    path_to_config = os.path.join(folder, "plot_config.json")
    with open(path_to_config) as f:
        config = json.load(f)

    title = request.GET.get('title', None)
    if title:
        config.update({"title": title})
    with open(path_to_config, 'w') as f:
        json.dump(config, f)
    call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "heatmap.R"), path_to_config]
    os.system(" ".join(call_list))
    return JsonResponse({"test":"ejejhe"})

def volcano_recalculate(request):
    id = request.GET.get('id', None)
    folder = os.path.join(MEDIA_ROOT, id)
    path_to_config = os.path.join(folder, "plot_config.json")
    with open(path_to_config) as f:
        config = json.load(f)

    title = request.GET.get('title', None)
    if title:
        config.update({"title": title})
    with open(path_to_config, 'w') as f:
        json.dump(config, f)
    call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "volcano.R"), path_to_config]
    os.system(" ".join(call_list))
    return JsonResponse({"test":"ejejhe"})


class Heatmap(FormView):
    #template_name = 'bench.html'
    #form_class = sRNABenchForm
    #success_url = reverse('photos:multi_start')

    def get_form_kwargs(self):
        '''This goes in the Update view'''
        kwargs = super(Heatmap, self).get_form_kwargs()  # put your view name in the super

        #kwargs["folder"] = self.request.path
        return kwargs

    def get(self, request,**kwargs):
        path = request.path
        folder = path.split("/")[-1]
        heatmap_url = os.path.join(MEDIA_URL,folder,"heatmap.html")



        return render(self.request, 'heatmap_template.html',
                      {"job_id": folder,
                       "heatmap_url": heatmap_url
                       })

class Volcano(FormView):
    #template_name = 'bench.html'
    #form_class = sRNABenchForm
    #success_url = reverse('photos:multi_start')

    def get_form_kwargs(self):
        '''This goes in the Update view'''
        kwargs = super(Volcano, self).get_form_kwargs()  # put your view name in the super

        #kwargs["folder"] = self.request.path
        return kwargs

    def get(self, request,**kwargs):
        path = request.path
        folder = path.split("/")[-1]
        volcano_url = os.path.join(MEDIA_URL,folder,"volcano.html")



        return render(self.request, 'volcano.html',
                      {"job_id": folder,
                       "volcano_url": volcano_url
                       })

