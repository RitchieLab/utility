import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import math
import numpy as np

###Formatting
gdat = pd.read_csv("~/Desktop/ACTG_GWAS_TEST.txt", sep=" ")
gdat['Tissue'] = "GWAS"
tdat = pd.read_csv("~/Desktop/ACTG_TWAS_TEST_map.txt", sep=" ")
dat = gdat.append(tdat, ignore_index=True)
dat[['CHR:BP', 'POS']] = dat['CHR:BP'].str.split(':', expand=True)
dat = dat[['PHE', 'SNP', 'CHR:BP', 'POS', 'A1:MAF', 'Interaction', 'N', 'pvalue', 'Tissue']]
dat.columns = ['PHE', 'SNP', 'CHR', 'POS', 'MAF', 'Interaction', 'N', 'pvalue', 'Tissue']
dat['CHR'] = dat.CHR.astype(int)
dat['POS'] = dat.POS.astype(int)
d_order = dat.sort_values(by=['CHR','POS'])
d_order['pos_index'] = range(1, len(d_order)+1)
#Get tick positions
maxRows = d_order.loc[d_order.groupby('CHR').pos_index.agg('idxmax')][['CHR', 'POS', 'pos_index']]
minRows = d_order.loc[d_order.groupby('CHR').pos_index.agg('idxmin')][['CHR', 'POS', 'pos_index']]
maxRows.columns = ['chrmax', 'posmax', 'posindmax']
minRows.columns = ['chrmin', 'posmin', 'posindmin']
lims = pd.concat([minRows.reset_index(), maxRows.reset_index()], axis=1)
lims['av'] = lims[['posindmax', 'posindmin']].mean(axis=1)
#Get every other row for background shape
limsshape = lims.iloc[1::2, :]
d_order['logp'] = -np.log(d_order['pvalue'])
ymin = math.floor(d_order['logp'].min())
ymax = math.ceil(d_order['logp'].max())
xmin = 0
xmax = d_order['pos_index'].max()

#Create shapes for chromosomes
sl = []
for i in range(int(len(limsshape.chrmin.unique()))):
    shapes = {'type': 'rect', 'x0': limsshape.iloc[i]['posindmin'], 'x1': limsshape.iloc[i]['posindmax'], 'y0': ymin-0.5, 'y1': ymax, 'fillcolor': '#d3d3d3', 'opacity': 0.3, 'layer': 'below', 'line':{'width': 0}}
    sl.append(shapes)

#Create color map for phenotypes
#colorsIdx = {'p_altastgr2w48': 'rgb(215,48,39)', 'p_altastgr3w48': 'rgb(215,148,39)', 'p_altastgr4w48': 'rgb(66, 134, 244)'}

#Set up the menus
app = dash.Dash('Hello World')
#Make these into a loop to get the label and values

#Get list of interaction terms for menu
il = []
for i in d_order.Interaction.unique():
   interactions = {'label': i, 'value': i}
   il.append(interactions)

#Get list of tissue terms for menu
tl = []
for i in d_order.Tissue.unique():
    if i != "GWAS":
        tissues = {'label': i, 'value': i}
        tl.append(tissues)


app.layout = html.Div([
    html.Div([
        dcc.Dropdown(
            id='interaction-dropdown',
            options = il
        ),
    ],style={'width': '48%', 'display': 'inline-block'}),

    html.Div([
        dcc.Dropdown(
            id='tissue-dropdown',
            options = tl,
        ),
    ],style={'width': '48%', 'display': 'inline-block'}),

    html.Div([
        dcc.Graph(id='twas-graph'),
        dcc.Graph(id='gwas-graph'),
    ], style={'width': '100%', 'display': 'inline-block'})
])

@app.callback(Output('gwas-graph', 'figure'), [Input('interaction-dropdown', 'value'), Input('tissue-dropdown', 'value')])
def update_graph(selected_dropdown_value, selected_dropdown_value2):
    #Here it would be possible to split df into 9 separate traces and return the list   
    df = d_order[(d_order.Interaction==selected_dropdown_value) & (d_order.Tissue=="GWAS")]
    #cols = df['PHE'].map(colorsIdx)
    dl = []
    for i in df.PHE.unique():
        dfsub = df[df.PHE==i]
        dff = {'x': dfsub.pos_index, 'y': dfsub.logp, 'mode': 'markers', 'text': dfsub.SNP + '\n' + dfsub.MAF, 'name': i}
        dl.append(dff)
    return {
        'data': dl,
        'layout': {'shapes': sl, 'title': 'GWAS Plot', 'showlegend': True, 'xaxis': {'range': [xmin, xmax], 'ticktext': lims.chrmin, 'tickvals': lims.av, 'showgrid': False}, 'yaxis': {'range': [ymin,ymax], 'title': '-log10(pvalue)'}, 'hovermode': 'closest'}
    }

@app.callback(Output('twas-graph', 'figure'), [Input('interaction-dropdown', 'value'), Input('tissue-dropdown', 'value')])
def update_graph(selected_dropdown_value, selected_dropdown_value2):
    #Here it would be possible to split df into 9 separate traces and return the list   
    df = d_order[(d_order.Interaction==selected_dropdown_value) & (d_order.Tissue==selected_dropdown_value2)]
    #cols = df['PHE'].map(colorsIdx)
    dl = []
    for i in df.PHE.unique():
        dfsub = df[df.PHE==i]
        dff = {'x': dfsub.pos_index, 'y': dfsub.logp, 'mode': 'markers', 'text': dfsub.SNP, 'name': i}
        dl.append(dff)
    return {
        'data': dl,
        'layout': {'shapes': sl, 'title': 'TWAS Plot', 'showlegend': True, 'xaxis': {'range': [xmin, xmax], 'ticktext': lims.chrmin, 'tickvals': lims.av, 'showgrid': False}, 'yaxis': {'range': [ymin,ymax], 'title': '-log10(pvalue)'}, 'hovermode': 'closest'}
    }

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})

if __name__ == '__main__':
    app.run_server()

