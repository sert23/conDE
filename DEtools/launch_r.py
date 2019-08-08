import os
import subprocess
import sys
import json


#Example line python3 launch_r.py config.json

#  exosome,exosome,exosome,exosome,exosome,exosome,cell,cell,cell,cell,cell,cell /shared/home_ugr/eap/ cucu6


with open(sys.argv[1], "r") as file:
    config = json.load(file)


scripts_dict = { "de_deseq.R" : "DESeq",
                 "de_deseq2.R" : "DESeq2",
                 "de_edger.R" : "edgeR",
                 "de_noiseq.R" : "NOISeq",
                 "wilcoxon-test.R" : "Wilcoxon",
                 "t-test.R" : "T-test",
                 "de_limma-voom.R": "limma-voom",
                 "de_limma-trend.R": "limma-trend"
                }
# print("This is the name of the script: ", sys.argv[0])
# print("Number of arguments: ", len(sys.argv))
# print("The arguments are: " , str(sys.argv))

scripts_folder = config["scripts_folder"]



rscripts = [ f for f in os.listdir(scripts_folder) if os.path.isfile(os.path.join(scripts_folder, f))]

outdir = config["outdir"]
if not os.path.exists(outdir):
    os.mkdir(outdir)

run_list = []
for script in rscripts:
    folder_path = os.path.join(outdir,scripts_dict.get(script))
    os.mkdir(folder_path)
    call_list = [ config["Rscript_path"],
                  os.path.join(scripts_folder, script),
                  config["input_matrix"],
                  config["matrixDesc"],
                  os.path.join(outdir,folder_path),
                  config["base_name"]
    ]
    print(" ".join(call_list))
    run_list.append(subprocess.Popen(call_list,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE))

exit_codes = [p.wait() for p in run_list]



