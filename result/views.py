from django.shortcuts import render
from django.views.generic import FormView, DetailView
from django.http import JsonResponse
from conDE.settings import MEDIA_ROOT, MEDIA_URL

import os

# Create your views here.

# def DEresult(request):
#     results = {}
#     if 'id' in request.GET:
#         job_id = request.GET['id']
#         results["job_id"] = job_id
#
#     return render(request, "DE_result.html", results)


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
        onlyfiles = []
        method_list = ["edgeR", "DESeq", "DESeq2","NOISeq"]
        plot_list = [["UpSet","UpSet"],["Barplot","Barplot"]]
        start_plot = "/Users/ernesto/PycharmProjects/conDE/DEtools/tests/test_dir/Barplot.png".replace(MEDIA_ROOT,MEDIA_URL)
        # if os.path.exists(os.path.join(MEDIA_ROOT,folder)):
        #     onlyfiles = [[f,os.path.join(MEDIA_URL, folder, f)] for f in listdir(os.path.join(MEDIA_ROOT,folder)) if
        #              os.path.isfile(os.path.join(os.path.join(MEDIA_ROOT, folder), f))]
        #
        #     for x in ["conf.txt","SRR_files.txt","URL_files.txt"]:
        #         if [x,os.path.join(MEDIA_URL, folder, x)] in onlyfiles:
        #             onlyfiles.remove([x,os.path.join(MEDIA_URL, folder, x)])
        #
        # else:
        #     onlyfiles = []
        #     os.mkdir(os.path.join(MEDIA_ROOT,folder))
        #     JobStatus.objects.create(job_name=folder+"_multi", pipeline_key=folder, job_status="not_launched",
        #                          start_time=datetime.datetime.now(),
        #                          all_files=" ",
        #                          modules_files=" ",
        #                          pipeline_type="multiupload",
        #                          )
        return render(self.request, 'DE_result.html',
                      {"job_id": folder,
                       "method_list" : method_list,
                       "plot_list" : plot_list,
                       "start_plot": start_plot
                       })

def test_graph_ajax(request):

    test_string = request.GET.get('test_string', None)
    test_directory = "/Users/ernesto/PycharmProjects/conDE/DEtools/tests/test_dir".replace(MEDIA_ROOT,MEDIA_URL)
    path = os.path.join(test_directory,test_string+".png")
    data = {
        'path': path
    }
    return JsonResponse(data)
