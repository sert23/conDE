from django.shortcuts import render, redirect
from conDE.settings import MEDIA_ROOT, MEDIA_URL, RSCRIPT_PATH, LOCAL_TEST, RPLOTS_PATH, SUB_SITE
import pandas as pd
import os
from django.http import JsonResponse
# Create your views here.
from django.views.generic import FormView, DetailView
import json
import random
import string
import shutil


def random_redirect(request):
    id = request.GET.get('id', None)
    plot = request.GET.get('plot', None)
    return(redirect(os.path.join(MEDIA_URL,id,plot)))

def make_random_url(id,plot):
    random_url = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    url = os.path.join(SUB_SITE,"plot","random",random_url)
    # return url+"?plot=" + plot + "&id=" + id
    return os.path.join(MEDIA_URL,id,plot+"?var="+random_url)

def rnd_str():
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))



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

    new_outdir = os.path.join(folder,"volcano",rnd_str())
    os.mkdir(new_outdir)
    with open(path_to_config) as f:
        config = json.load(f)

    title = request.GET.get('title', None)
    if title:
        config.update({"title": title})
    pval = request.GET.get('pval', None)
    if title:
        config.update({"pval": pval})
    FC = request.GET.get('FC', None)
    if title:
        config.update({"FC": FC})
    config.update({"folder":new_outdir})
    with open(path_to_config, 'w') as f:
        json.dump(config, f)
    call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "volcano.R"), path_to_config]
    os.system(" ".join(call_list))
    return JsonResponse( {"new_url" : os.path.join(new_outdir.replace(MEDIA_ROOT,MEDIA_URL),"volcano.html")
                                                                                                        })




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
        #get id from link
        folder = path.split("/")[-1]
        folder_path = os.path.join(MEDIA_ROOT,folder)
        init_json = os.path.join(folder_path,"init_plot_config.json")
        current_json = os.path.join(folder_path,"plot_config.json")
        #copy initial json to "restart" the plot
        shutil.copy(init_json,current_json)
        #clear folder of older jobs
        shutil.rmtree(os.path.join(folder_path,"volcano"))
        os.mkdir(os.path.join(folder_path,"volcano"))
        #launch plotting job
        call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "volcano.R"), current_json]
        os.system(" ".join(call_list))
        volcano_url = make_random_url(folder,"volcano.html")

        return render(self.request, 'volcano.html',
                      {"job_id": folder,
                       "volcano_url": volcano_url
                       })

