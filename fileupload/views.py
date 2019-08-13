from django.shortcuts import render, redirect
import os
from django.views.generic import FormView
from conDE.settings import MEDIA_ROOT, MEDIA_URL, RSCRIPT_PATH, LOCAL_TEST, RPLOTS_PATH, SUB_SITE,RSCRIPTS_FOLDER,PYTHON_PATH,PATH_TO_DE_LAUNCHER
import random
import string
from django.http import JsonResponse
from .forms import PhotoForm
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.core.files.storage import FileSystemStorage
import plotly.graph_objs as go
from plotly.offline import plot
import pandas as pd
import csv
from django.urls import reverse_lazy
from shutil import copy
import json



# Create your views here.

def col_num(input_file,separator):

    with open(input_file) as f:
        reader = csv.reader(f, delimiter=separator, skipinitialspace=True)
        first_row = next(reader)

        return first_row

def generate_folder():
    stay = True
    while stay:
        random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        folder = os.path.join(MEDIA_ROOT,random_string)
        if not os.path.exists(folder):
            os.mkdir(folder)
            stay = False
            return random_string

def removeColumns(input_matrix,to_drop):
    input_table = pd.read_csv(input_matrix, sep="\t")
    input_table = input_table.drop(to_drop,axis=1)
    input_table.to_csv(input_matrix, sep="\t", index=False)

def example_table():
    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'

    fig = go.Figure(data=[go.Table(
        header=dict(
            values=['<b>EXPENSES</b>', '<b>Q1</b>', '<b>Q2</b>', '<b>Q3</b>', '<b>Q4</b>'],
            line_color='darkslategray',
            fill_color=headerColor,
            align=['left', 'center'],
            font=dict(color='white', size=12)
        ),
        cells=dict(
            values=[
                ['Salaries', 'Office', 'Merchandise', 'Legal', '<b>TOTAL</b>'],
                [1200000, 20000, 80000, 2000, 12120000],
                [1300000, 20000, 70000, 2000, 130902000],
                [1300000, 20000, 120000, 2000, 131222000],
                [1400000, 20000, 90000, 2000, 14102000]],
            line_color='darkslategray',
            # 2-D list of colors for alternating rows
            fill_color=[[rowOddColor, rowEvenColor, rowOddColor, rowEvenColor, rowOddColor] * 5],
            align=['left', 'center'],
            font=dict(color='darkslategray', size=11)
        ))
    ])
    # plot(fig, filename="/Users/ernesto/PycharmProjects/conDE/upload/AA99/yu.html",  show_link=False, auto_open=False, include_plotlyjs=True)
    return plot(fig,  show_link=False, auto_open=False, include_plotlyjs=True, output_type='div')

def preview_table(job_id, separator="\t", has_header=True, first_missing=False):
    target_file = os.path.join(MEDIA_ROOT, job_id, "uploaded", "input.matrix")
    if has_header:
        if first_missing:
            headers = ["name"] + col_num(target_file, separator)
            input_table = pd.read_csv(target_file, sep=separator, names=headers)
        else:
            input_table = pd.read_csv(target_file,sep=separator)
    else:
        cols = col_num(target_file, separator)
        col_n = len(cols)
        a = range(col_n-1)
        b = (col_n-1)*["sample-"]
        headers = ["name"] + ["{}{:02}".format(b_, a_) for a_, b_ in zip(a, b)]
        input_table = pd.read_csv(target_file, sep=separator, names=headers)

    # pre_table = input_table.head(5).iloc[:, : 6]
    pre_table = input_table.head(5)
    pre_table.columns.values[0] = "name"
    # print(pre_table)
    output_file = os.path.join(MEDIA_ROOT, job_id, "uploaded", "parsed.matrix")
    input_table.to_csv(output_file,sep="\t", index=False)
    headerColor = 'grey'
    rowEvenColor = 'lightgrey'
    rowOddColor = 'white'
    fig = go.Figure(data=[go.Table(
        header=dict(
            values= pre_table.columns ,
            line_color='darkslategray',
            fill_color=headerColor,
            align=['left', 'center'],
            font=dict(color='white', size=12)
        ),
        cells=dict(
            values= pre_table.T.values.tolist(),
            line_color='darkslategray',
            # 2-D list of colors for alternating rows
            # fill_color=[[rowOddColor, rowEvenColor, rowOddColor, rowEvenColor, rowOddColor] * 5],
            align=['left', 'center'],
            font=dict(color='darkslategray', size=11)
        ))
    ],
        layout= go.Layout(
            autosize=False,
            width=10 * pre_table.name.map(len).max()*len(pre_table.columns))
            #10*pre_table.columns.map(len).max()*(len(pre_table.columns)-1)) for columns width
    )

    return plot(fig,  show_link=False, auto_open=False, include_plotlyjs=False, output_type='div')
    # return plot(fig,  show_link=False, auto_open=False, include_plotlyjs=True, output_type='div')



class Upload(FormView):
    #template_name = 'bench.html'
    #form_class = sRNABenchForm
    #success_url = reverse('photos:multi_start')

    # @method_decorator(csrf_exempt)
    def dispatch(self, request, *args, **kwargs):
        return super(Upload, self).dispatch(request, *args, **kwargs)

    def get_form_kwargs(self):
        '''This goes in the Update view'''
        kwargs = super(Upload, self).get_form_kwargs()  # put your view name in the super

        #kwargs["folder"] = self.request.path
        return kwargs

    def get(self, request,**kwargs):
        path = request.path
        folder = generate_folder()
        os.mkdir(os.path.join(MEDIA_ROOT, folder,"uploaded"))
        header_list = ["sample1","sample2","sample3","sample4","sample5"]
        button_link = reverse_lazy("groups") + "/" +folder
        # return render(self.request, 'tagsinput_start.html', {"request_path":path, "job_id" : folder,
        #                                                      "form" : PhotoForm, "header_list":header_list })
        return render(self.request, 'simple_upload.html', {"request_path":path, "job_id" : folder, "form" : PhotoForm,
                                                           "header_list": header_list, "button_link": button_link})
        # return render(self.request, 'multiupload.html', {'file_list': onlyfiles, "request_path":path, "form": MultiURLForm })
        #return render(self.request, 'multiupload.html', {'file_list': [os.path.join(MEDIA_ROOT,folder),os.path.join(MEDIA_ROOT,folder)]})

    # @method_decorator(csrf_exempt)
    def post(self, request):
        sep_dict = {"tab":"\t", "comma":",",  "semicolon":";" }
        job_id = request.POST.get("job_id")
        print(job_id)
        header = bool(request.POST.get("header"))
        separator = request.POST.get("separator")
        first_missing = bool(request.POST.get("missing"))
        sep = sep_dict.get(separator)
        target_file = os.path.join(MEDIA_ROOT,job_id,"uploaded","input.matrix")
        if os.path.exists(target_file):
            os.remove(target_file)
        group1 = request.POST.get("id_group1").replace(" ","_")
        group2 = request.POST.get("id_group2").replace(" ","_")
        if group1 and group2:
            with open(os.path.join(MEDIA_ROOT,job_id,"groups"),"w") as gf:
                gf.write(group1+"\n")
                gf.write(group2+"\n")
        if "file" in self.request.FILES:
            for f in self.request.FILES.getlist('file'):
                fs = FileSystemStorage(location=MEDIA_ROOT)
                filename = fs.save(target_file, f)
        table = preview_table(job_id,sep, header,first_missing)
        return JsonResponse({'error': False, 'message': 'Uploaded Successfully', "table":table})

                # return JsonResponse()
                # return JsonResponse(data)

        # else:
        #     dfolder = os.path.join(MEDIA_ROOT, folder)
        #     form = MultiURLForm(self.request.POST, self.request.FILES, dest_folder= dfolder)
        #     form.is_valid()
        #     form.clean()
        #     if form.is_valid():
        #         url = reverse('multi:multi_launch') + folder
        #     #url = reverse('photos:multi_start')
        #         return redirect(url)


class ChooseGroups(FormView):
    #template_name = 'bench.html'
    #form_class = sRNABenchForm
    #success_url = reverse('photos:multi_start')

    # @method_decorator(csrf_exempt)
    def dispatch(self, request, *args, **kwargs):
        return super(ChooseGroups, self).dispatch(request, *args, **kwargs)

    def get_form_kwargs(self):
        '''This goes in the Update view'''
        kwargs = super(ChooseGroups, self).get_form_kwargs()  # put your view name in the super

        #kwargs["folder"] = self.request.path
        return kwargs

    def get(self, request,**kwargs):
        path = request.path
        folder = path.split("/")[-1]
        initial_file = os.path.join(MEDIA_ROOT, folder, "uploaded","parsed.matrix")
        target_file = os.path.join(MEDIA_ROOT, folder, "input.matrix")
        copy(initial_file, target_file)
        header_list = col_num(target_file, "\t")[1:]
        with open(os.path.join(MEDIA_ROOT, folder, "groups"), "r") as gf:
            lines = gf.readlines()
            group1 = lines[0].rstrip().replace("_", " ")
            group2 = lines[1].rstrip().replace("_", " ")


        # print("hello")
        return render(self.request, 'choose_groups.html', {"request_path":path, "job_id" : folder,
                                                           "header_list": header_list, "group1":group1, "group2":group2
                                                           })

    def post(self, request):
        path = request.path
        folder = path.split("/")[-1]
        target_file = os.path.join(MEDIA_ROOT, folder, "input.matrix")
        header_list = col_num(target_file, "\t")[1:]
        with open(os.path.join(MEDIA_ROOT, folder, "groups"), "r") as gf:
            lines = gf.readlines()
            group1_str = lines[0].rstrip()
            group2_str = lines[1].rstrip()

        job_id = request.POST.get("job_id")
        group1 = request.POST.get("id_group1").split(",")
        group2 = request.POST.get("id_group2").split(",")
        group1_index = [i for i, item in enumerate(header_list) if item in group1]
        group2_index = [i for i, item in enumerate(header_list) if item in group2]
        if len(header_list) > len(group2_index+group1_index):
            to_drop = [x for x in header_list if x not in group1+group2]
            removeColumns(target_file,to_drop)
            header_list = col_num(target_file, "\t")[1:]

        gr_list = [group1_str if x in group1 else x for x in header_list]
        gr_list = [group2_str if x in group2 else x for x in gr_list]

        config_dict = {"outdir": os.path.join(MEDIA_ROOT,folder,"de"),
                       "input_matrix" : target_file,
                       "Rscript_path" : RSCRIPT_PATH,
                       "base_name": "_",
                       "matrixDesc": ",".join(gr_list),
                       "scripts_folder": RSCRIPTS_FOLDER}
        config_path = os.path.join(MEDIA_ROOT,folder,"de_config.json")

        with open(config_path,"w") as cf:
            json.dump(config_dict, cf)

        launch_string = " ".join([PYTHON_PATH,PATH_TO_DE_LAUNCHER, config_path])
        with open(os.path.join(MEDIA_ROOT,folder,"test.out"), "w") as text_file:
            text_file.write(launch_string)
        os.system(launch_string)

        redirection = reverse_lazy("result") + "/" +folder

        return redirect(redirection)

