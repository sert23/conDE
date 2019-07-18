from django.shortcuts import render
from django.views.generic import FormView, DetailView
from django.http import JsonResponse
from conDE.settings import MEDIA_ROOT, MEDIA_URL, RSCRIPT_PATH, LOCAL_TEST, RPLOTS_PATH
import pandas as pd
import os
import subprocess
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot
import json
from random import randrange


# Create your views here.

# def DEresult(request):
#     results = {}
#     if 'id' in request.GET:
#         job_id = request.GET['id']
#         results["job_id"] = job_id
#
#     return render(request, "DE_result.html", results)

def update_json(folder_path, FC=None,pval=None,iset=None,methods=None):
    path_to_config = os.path.join(folder_path,"config.json")
    with open(path_to_config) as f:
        config = json.load(f)
    changed = False
    if FC != config["FC"]:
        config["FC"] = FC
        changed = True
    if pval != config["pval"]:
        config["pval"] = pval
        changed = True

    config["set"] = iset

    if set(methods) != set(config["methods"]):
        config["methods"] = methods
        changed = True

    with open(path_to_config, 'w') as f:
        json.dump(config, f)
    return changed

def select_DE(path,pval,FC,out):
    with open(path,"r") as in_path:
        df = pd.read_csv(in_path,sep="\t")

    FC = float(FC)
    pval = float(pval)
    if FC < 1:
        FC= 1/float(FC)
    k1 = df.loc[(df.pvalue < pval) & ((df.FoldChange > FC) | (df.FoldChange < 1/float(FC)))]
    k1.to_csv(out, sep='\t')
    return k1['name'].tolist()

def calculate_consensus(directory,methods,pval,FC):
    import itertools
    import pandas
    con_list=[]
    for method in methods:
        out_dir = os.path.join(directory,"consensus",method+".txt")
        de_dir = os.path.join(directory,"de",method)
        in_file = [f for f in os.listdir(de_dir) if f.endswith("allGenes.csv")][0]
        in_path = os.path.join(directory,"de",method,in_file)
        con_list.append(select_DE(in_path,pval,FC,out_dir))
    merged = list(itertools.chain.from_iterable(con_list))
    ps = pandas.Series(merged)
    counts = ps.value_counts()
    with open(os.path.join(directory,"consensus.txt"),"w") as con:
        for i, val in counts.iteritems():
            if val >= len(methods):
                con.write(i+"\n")

def draw_barplot(path_to_config):
    with open(path_to_config) as f:
        config = json.load(f)
    n_over = {}
    n_under = {}
    for method in config["methods"]:
        in_file = os.path.join(config["folder"],"consensus", method+".txt")
        with open(in_file, "r") as in_path:
            df = pd.read_csv(in_path, sep="\t")
        n_under[method] = len(df.loc[(df.FoldChange < 1 )]["name"].tolist())
        n_over[method] = len(df.loc[(df.FoldChange > 1 )]["name"].tolist())

    data = []
    if (config["set"] =="All") or (config["set"] =="Over"):
        over = go.Bar(
            x= list(n_over.keys()) ,
            y= list(n_over.values()),
            marker_color= '#2ca02c',
            marker=dict(line=dict(width = 1 ,color="black")),

            name='Number of overexpressed')
        data.append(over)
    if (config["set"] =="All") or (config["set"] =="Under"):
        under = go.Bar(
            x= list(n_under.keys()) ,
            y= list(n_under.values()),
            marker_color='crimson',
            marker=dict(line=dict(width=1, color="black")),
            name='Number of underexpressed')
        data.append(under)
    layout = go.Layout(
            barmode='stack'
        )
    fig = go.Figure(data=data, layout=layout)
    out_path = os.path.join(config["folder"],"plots","BarPlot.html")
    plot(fig, filename=out_path,show_link=False, auto_open=False, include_plotlyjs=True)
    print(out_path)
    # return out_path
    # return plot(fig, show_link=False, auto_open=False, output_type='div', include_plotlyjs=True)



def draw_plots(path_to_config):
    wait_list = []
    call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH,"upset.R"),  path_to_config]
    draw_barplot(path_to_config)
    if not LOCAL_TEST:
        wait_list.append(subprocess.Popen(call_list,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE))
        call_list = [RSCRIPT_PATH, os.path.join(RPLOTS_PATH,"VennD.R"), path_to_config]
        wait_list.append(subprocess.Popen(call_list,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE))
        exit_codes = [p.wait() for p in wait_list]
    folders = path_to_config.split("/")
    temp = "/".join(folders[:-1])
    with open(os.path.join(temp,"temp"),"w") as tmp:
        tmp.write(" ".join(call_list))



def get_plot_content(request):
    data = {}
    plot = request.GET.get('plot', None).replace(" ","_")
    id = request.GET.get('id', None)
    plot_path = os.path.join(MEDIA_ROOT,id,"plots",plot+".jpg")
    if os.path.exists(plot_path):
        media_plot = plot_path.replace(MEDIA_ROOT,MEDIA_URL)
        div_content = ' <div class="col-lg-12"> <img src="' + media_plot + '?'+ str(randrange(500)) +'" style="width:100%;height:100%;padding:1px;border:thin solid black;">  </div> '
    else:
        div_content = '<iframe style="border-style:solid;" src='+ os.path.join(MEDIA_ROOT,id,"plots","BarPlot.html").replace(MEDIA_ROOT,MEDIA_URL)  +' width="100%" height="500" allowfullscreen></iframe>'

    data["div_content"] = div_content
    return JsonResponse(data)

class DEresult(FormView):
    #template_name = 'bench.html'
    #form_class = sRNABenchForm
    #success_url = reverse('photos:multi_start')

    def get_form_kwargs(self):
        '''This goes in the Update view'''
        kwargs = super(DEresult, self).get_form_kwargs()  # put your view name in the super

        #kwargs["folder"] = self.request.path
        return kwargs

    def get(self, request,**kwargs):
        path = request.path
        folder = path.split("/")[-1]
        folder_path = os.path.join(MEDIA_ROOT,folder)
        de_path = os.path.join(folder_path,"de")
        method_list = [f for f in os.listdir(de_path) if os.path.isdir(os.path.join(de_path, f))]
        method_list.sort()
        update_json(folder_path, FC=2, pval=0.05, iset="All", methods=method_list)
        calculate_consensus(folder_path, method_list, 0.05, 2)
        # method_list = ["edgeR", "DESeq", "DESeq2","NOISeq"]
        plot_list = [["UpSet","UpSet"],["Barplot","Barplot"],["Venn","Venn Diagram"]]
        set_list = [["All","All DE genes"],["Over","Overexpressed only"],["Under","Underexpressed only"]]
        calculate_consensus(folder_path,method_list,0.05,2)
        draw_plots(os.path.join(folder_path,"config.json"))
        start_plot = os.path.join(folder_path,"plots","UpSet.jpg").replace(MEDIA_ROOT,MEDIA_URL)


        return render(self.request, 'DE_result.html',
                      {"job_id": folder,
                       "method_list" : method_list,
                       "plot_list" : plot_list,
                       "start_plot": start_plot,
                       "set_list": set_list,
                       })



def ajax_recalculate(request):
    mock = {}
    id = request.GET.get('id', None)
    FC = request.GET.get('FC', None)
    pval = request.GET.get('pval', None)
    methods = request.GET.get('methods', None)
    methods = methods.split(",")
    iset = request.GET.get('set', None)
    folder = os.path.join(MEDIA_ROOT, id )
    if update_json(folder, FC=FC,pval=pval,iset=iset,methods=methods):
        calculate_consensus(folder, methods, pval, FC)
    draw_plots(os.path.join(folder, "config.json"))
    return JsonResponse(mock)
