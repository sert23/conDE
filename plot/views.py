from django.shortcuts import render, redirect
from conDE.settings import MEDIA_ROOT, MEDIA_URL, RSCRIPT_PATH, LOCAL_TEST, RPLOTS_PATH, SUB_SITE,PYTHON_PATH
import pandas as pd
import os
from django.http import JsonResponse
# Create your views here.
from django.views.generic import FormView, DetailView
import json
import random
import string
import shutil
from django.urls import reverse_lazy



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

def restrict_to_consensus(input_matrix,input_consensus):
    df = pd.read_csv(input_matrix, sep="\t")
    consensus_list = []
    with open(input_consensus, "r") as consensus_file:
        lines = consensus_file.readlines()
    for line in lines:
        consensus_list.append(line.rstrip())
    defin = df[df['name'].isin(consensus_list)]
    defin = defin.drop(["pvalue","padj"],axis=1)
    with open(input_matrix,"w") as out_f:
        defin.to_csv(input_matrix, sep='\t', index=False)



def new_plot_view(request):
    old_id = request.GET.get('job_id', None)
    plot = request.GET.get('id_plot', None)
    method = request.GET.get('id_plot_method', None)
    old = True
    while old:
        new_id = rnd_str()
        new_path = os.path.join(MEDIA_ROOT,new_id)
        if not os.path.exists(new_path):
            os.mkdir(new_path)
            old = False
    orig_tab = os.path.join(MEDIA_ROOT,old_id,"de",method,"allGenes.csv")
    dest_tab = os.path.join(new_path,"input.matrix")
    shutil.copy(orig_tab,dest_tab)
    if request.GET.get('consensus', None):
        restrict_to_consensus(dest_tab, os.path.join(MEDIA_ROOT,old_id,"consensus.txt"))

    with open(os.path.join(MEDIA_ROOT,old_id,"de_config.json")) as f:
        original_config = json.load(f)
    original_config.update({"input_matrix":dest_tab,
                             "title": " ", "top_n": "20",
                            "folder": os.path.join(MEDIA_ROOT,new_id)
                            })
    with open(os.path.join(MEDIA_ROOT, new_id, "init_plot_config.json"), 'w') as f:
        json.dump(original_config, f)

    if plot == "heatmap":
        return redirect(reverse_lazy("heatmap") + new_id)
    if plot == "volcano":
        return redirect(reverse_lazy("heatmap") + new_id)
    if plot == "PCA":
        return redirect(reverse_lazy("pca") + new_id)


    print(old_id)
    print(plot)
    print(new_id)
    print(method)
    #calculate random string
    #make initial config
    #copy expression matrix



def heatmap_recalculate(request):
    id = request.GET.get('id', None)
    folder = os.path.join(MEDIA_ROOT, id)
    path_to_config = os.path.join(folder, "plot_config.json")
    shutil.rmtree(os.path.join(folder, "heatmap"))
    os.mkdir(os.path.join(folder, "heatmap"))
    new_outdir = os.path.join(folder, "heatmap", rnd_str())
    os.mkdir(new_outdir)

    with open(path_to_config) as f:
        config = json.load(f)

    group1_name = request.GET.get('group1_name', None)
    if group1_name:
        print("")
    group2_name = request.GET.get('group2_name', None)
    if group2_name:
        print("")
    ntop = request.GET.get('ntop', None)
    if ntop:
        config.update({'ntop': ntop})
    sortBy = request.GET.get('sortBy', None)
    if sortBy:
        config.update({'sortBy': sortBy})
    sortSense = request.GET.get('sortSense', None)
    if sortSense:
        config.update({'sortSense': sortSense})

    iset = request.GET.get('set', None)
    if iset:
        config.update({"set": iset})
    title = request.GET.get('title', None)
    if title:
        config.update({"title": title})
    pval = request.GET.get('pval', None)
    if title:
        config.update({"pval": pval})
    FC = request.GET.get('FC', None)
    if title:
        config.update({"FC": FC})
    config.update({"folder": new_outdir})
    with open(path_to_config, 'w') as f:
        json.dump(config, f)
    call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "heatmap.R"), path_to_config]
    os.system(" ".join(call_list))
    return JsonResponse({"new_url": os.path.join(new_outdir.replace(MEDIA_ROOT, MEDIA_URL), "heatmap.html")
                         })


def volcano_recalculate(request):
    id = request.GET.get('id', None)
    folder = os.path.join(MEDIA_ROOT, id)
    path_to_config = os.path.join(folder, "plot_config.json")

    shutil.rmtree(os.path.join(folder, "volcano"))
    os.mkdir(os.path.join(folder, "volcano"))

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

def header_finder(input_matrix, input_list):
    with open(input_matrix,"r") as f:
        first_line = f.readline().rstrip()
        headers = first_line.split("\t")
        return [x for x in input_list if x in headers]


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
        folder_path = os.path.join(MEDIA_ROOT,folder)
        init_json = os.path.join(folder_path, "init_plot_config.json")
        current_json = os.path.join(folder_path, "plot_config.json")
        # copy initial json to "restart" the plot
        shutil.copy(init_json, current_json)
        group1=""
        group2=""

        if os.path.exists(os.path.join(folder_path, "groups")):
            with open(os.path.join(folder_path,"groups"),"r") as gf:
                lines = gf.readlines()
                group1 = lines[0].rstrip().replace("_"," ")
                group2 = lines[1].rstrip().replace("_"," ")
        addFC = False
        addpval = False
        headers_found = header_finder(os.path.join(folder_path,"input.matrix"),["FoldChange","pvalue"])
        if "FoldChange" in headers_found:
            addFC=True
        if "pvalue" in headers_found:
            addpval=True
        if os.path.exists(os.path.join(folder_path, "heatmap")):
            shutil.rmtree(os.path.join(folder_path, "heatmap"))
        os.mkdir(os.path.join(folder_path, "heatmap"))
        if not os.path.exists(os.path.join(MEDIA_ROOT,folder,"heatmap.html")):
            call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "heatmap.R"), current_json]
            os.system(" ".join(call_list))

        heatmap_url = os.path.join(MEDIA_URL,folder,"heatmap.html")

        return render(self.request, 'heatmap_template.html',
                      {"job_id": folder,
                       "heatmap_url": heatmap_url,
                       "group1_name": group1,
                       "group2_name": group2,
                       "addFC": addFC,
                       "addpval": addpval,
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

class Pca(FormView):
    #template_name = 'bench.html'
    #form_class = sRNABenchForm
    #success_url = reverse('photos:multi_start')

    def get_form_kwargs(self):
        '''This goes in the Update view'''
        kwargs = super(Pca, self).get_form_kwargs()  # put your view name in the super

        #kwargs["folder"] = self.request.path
        return kwargs

    def get(self, request,**kwargs):
        path = request.path
        folder = path.split("/")[-1]
        folder_path = os.path.join(MEDIA_ROOT,folder)
        init_json = os.path.join(folder_path, "init_plot_config.json")
        current_json = os.path.join(folder_path, "plot_config.json")
        # copy initial json to "restart" the plot
        shutil.copy(init_json, current_json)
        group1=""
        group2=""

        if os.path.exists(os.path.join(folder_path, "groups")):
            with open(os.path.join(folder_path,"groups"),"r") as gf:
                lines = gf.readlines()
                group1 = lines[0].rstrip().replace("_"," ")
                group2 = lines[1].rstrip().replace("_"," ")
        addFC = False
        addpval = False
        headers_found = header_finder(os.path.join(folder_path,"input.matrix"),["FoldChange","pvalue"])
        if "FoldChange" in headers_found:
            addFC=True
        if "pvalue" in headers_found:
            addpval=True
        if os.path.exists(os.path.join(folder_path, "PCA")):
            shutil.rmtree(os.path.join(folder_path, "PCA"))
        os.mkdir(os.path.join(folder_path, "PCA"))
        if not os.path.exists(os.path.join(MEDIA_ROOT,folder,"PCA.html")):
            # call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH, "heatmap.R"), current_json]
            call_list = [PYTHON_PATH, os.path.join(RPLOTS_PATH, "PCA.py"), current_json]
            os.system(" ".join(call_list))

        plot_url = os.path.join(MEDIA_URL,folder,"PCA.html")

        return render(self.request, 'PCA_template.html',
                      {"job_id": folder,
                       "plot_url": plot_url,
                       "group1_name": group1,
                       "group2_name": group2,
                       "addFC": addFC,
                       "addpval": addpval,
                       })
