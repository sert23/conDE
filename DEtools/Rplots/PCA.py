import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from plotly.offline import plot
import numpy
import plotly.graph_objs as go
import json
import sys
import os

def plot_PCA(input_df, PCA_obj,config_dict):
    data = []
    uniq_groups = pd.unique(input_df.groups.values)

    for g in uniq_groups:
        k1 = input_df.loc[(input_df.groups == g)]
        trace = dict(
        type='scatter',
        x=k1.PC1.values,
        y=k1.PC2.values,
        mode='markers',
        name=g,
        marker=dict(
            # color=col,
            size=12,
            line=dict(
                color='rgba(217, 217, 217, 0.14)',
                width=0.5),
            opacity=0.8)
    )
        data.append(trace)

    # layout = go.Layout(
    #     margin=go.layout.Margin(
    #         r=500,
    #     ),
    #
    #     title="")

    layout = go.Layout(
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=100,
            t=100,
            pad=4
        ),
        title=config_dict.get("title"," "),
        font=dict(size=18),
        # autosize=False,
        # height=650,
        # width=1150,
        xaxis=dict(
            automargin=True,
            title='PC1',
            # tick0=0,
            # dtick=2,
        ),
        yaxis=dict(
            # type='log',
            automargin=True,
            # ticksuffix='%',
            # tickprefix="   ",
            title='PC2')
    )



    fig = dict(data=data, layout=layout)
    # div = plot(fig, output_type="div", show_link=False, auto_open=False, include_plotlyjs=True)
    plot(fig, filename=os.path.join(config_dict.get("folder"),"PCA.html"), show_link=False, auto_open=False, include_plotlyjs=True)
    # return div

input_json = sys.argv[1]
# input_json = "/Users/ernesto/PycharmProjects/conDE/upload/HRNYYDATFZ/plot_config.json"

with open(input_json, "r") as file:
    config = json.load(file)

groups = config["matrixDesc"]
input_file = config["input_matrix"]

df = pd.read_csv(
    filepath_or_buffer=input_file,
    sep='\t')

pval = config.get("pval")
if pval:
    pval = float(pval)
FC = config.get("FC")
if FC:
    FC = float(FC)
else:
    FC = float(1)
cset = config.get("set")
ntop = config.get("top_n")
if ntop:
    ntop = int(ntop)

if pval and "pval" in df.columns:
    df = df.loc[(df.pvalue < pval) & ((df.FoldChange > FC) | (df.FoldChange < 1/float(FC)))]
else:
    df = df.loc[((df.FoldChange > FC) | (df.FoldChange < 1 / float(FC)))]

if cset == "Over":
    df = df.loc[df.FoldChange > 1]
elif cset == "Under":
    df = df.loc[df.FoldChange < 1 ]

# df = df.head(ntop)

df = df.drop(["FoldChange","log2FoldChange","pvalue","padj"],axis=1, errors="ignore")
df_t = df.T
df_t.columns = df_t.iloc[0]
df_t = df_t.drop(df_t.index[0])

X = df_t.iloc[:,1:].values #features

y = groups.split(",") #class
pca = PCA()
xt = pca.fit_transform(X)
to_plot = [item[[0,1]] for item in xt]
print(dir(pca))
# print(pca.explained_variance_ratio_)

print(numpy.array(to_plot)[:,0])
print(numpy.array(to_plot)[:,1])
labels = list(df_t.index)
print(labels)
print("here")
print(y)
d = {'PC1': numpy.array(to_plot)[:,0] , 'PC2': numpy.array(to_plot)[:,1], 'groups':y, "labels": labels}
idf = pd.DataFrame(data=d)
plot_PCA(idf,pca, config)

