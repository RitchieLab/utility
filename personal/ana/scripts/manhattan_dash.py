import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import math
import numpy as np

###Formatting
gdat = pd.read_csv("~/Desktop/projects/plotly_dash_ACTG/ACTG_GWAS_TEST.txt", sep=" ")
gdat[['A1', 'MAF']] = gdat['A1:MAF'].str.split(':', expand=True)
gdat['Tissue'] = "GWAS"
gdat['TissueCategory']= "GWAS"
gdat['MAF'] = pd.to_numeric(gdat['MAF'])
gdat['MAF2'] = np.where(gdat['MAF']>0.5, 1-gdat['MAF'], gdat['MAF'])
gdat = gdat[['PHE', 'SNP', 'CHR:BP', 'A1', 'MAF2', 'Interaction', 'N', 'pvalue', 'Tissue', 'TissueCategory']]
gdat.columns = ['PHE', 'SNP', 'CHR:BP', 'A1', 'MAF', 'Interaction', 'N', 'pvalue', 'Tissue', 'TissueCategory']
tdat = pd.read_csv("~/Desktop/projects/plotly_dash_ACTG/ACTG_TWAS_TEST_map2.txt", sep=" ")
dat = gdat.append(tdat, ignore_index=True)
dat[['CHR', 'POS']] = dat['CHR:BP'].str.split(':', expand=True)
dat = dat[['PHE', 'SNP', 'CHR', 'POS', 'A1', 'MAF', 'Interaction', 'N', 'pvalue', 'Tissue', 'TissueCategory']]
dat['CHR'] = dat.CHR.astype(int)
dat['POS'] = dat.POS.astype(int)
dat['MAF'] = dat.MAF.astype(float)
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

#Hover
#d_order['Hover'] = np.where(d_order['Tissue']=='GWAS', d_order['SNP'] + "\nN:" + d_order['N'].map(str) + "\nMAF:" + d_order['MAF'].map(str) + "\nGene:" + d_order['Gene'].map(str), "N:" + d_order['N'].map(str) + "\nTissue:" + d_order['Tissue'] + "\nGene:" + d_order['Gene'].map(str))

#Create shapes for chromosomes
sl = []
for i in range(int(len(limsshape.chrmin.unique()))):
    shapes = {'type': 'rect', 'x0': limsshape.iloc[i]['posindmin'], 'x1': limsshape.iloc[i]['posindmax'], 'y0': ymin-0.5, 'y1': ymax, 'fillcolor': '#d3d3d3', 'opacity': 0.3, 'layer': 'below', 'line':{'width': 0}}
    sl.append(shapes)

#Create color map for phenotypes
#colorsIdx = {'p_altastgr2w48': 'rgb(215,48,39)', 'p_altastgr3w48': 'rgb(215,148,39)', 'p_altastgr4w48': 'rgb(66, 134, 244)'}

#Set up the menus
app = dash.Dash(
    meta_tags=[
            {"name": "viewport", "content": "width=device-width, initial-scale=1"}
    ]
)
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

#Create dictionary for tissues and tissue categories
tcdict = d_order.groupby('TissueCategory')['Tissue'].unique().apply(list).to_dict()
tcdict.pop('GWAS', None)

#Get list of tissue categories for menu
cl = []
for i in d_order.TissueCategory.unique():
    if i != "GWAS":
        categories = {'label' : i, 'value': i}
        cl.append(tissues)

#Get marker positions for slider based on max MAF
if math.ceil(d_order['MAF'].max()*10)/10 > 0.5:
   markers={0: str(0), 0.1: str(0.1), 0.2: str(0.2), 0.3: str(0.3), 0.4: str(0.4), 0.5: str(0.5), 0.6: str(0.6), 0.7: str(0.7), 0.8: str(0.8), 0.9: str(0.9), 1: str(1)}
   maxmaf = 1
else:
   markers={0: str(0), 0.05: str(0.05), 0.1: str(0.1), 0.15: str(0.15), 0.2: str(0.2), 0.25: str(0.25), 0.3: str(0.3), 0.35: str(0.35), 0.4: str(0.4), 0.45: str(0.45), 0.5: str(0.5)}
   maxmaf = 0.5

app.layout = html.Div([
    html.Div([
        dcc.Dropdown(
            id='interaction-dropdown',
            options = il,
            placeholder = "Select an interaction term..."
        ),
    ],style={'width': '49%', 'display': 'inline-block'}),
    #html.H1(children='MAF', style={'width': '45%', 'display': 'inline-block', 'font-size': '18px', 'textAlign': 'center'}), 

    html.Div([
        dcc.RangeSlider(
            id='maf-slider',
            min=0,
            max=maxmaf,
            step=maxmaf/10,
            marks=markers,
            value=[0,maxmaf],
            allowCross=False
        )
    ],style={'width': '49%', 'display': 'inline-block', 'verticalAlign': 'top', 'horizontalAlign': 'right'}),

    html.Div([
        dcc.Dropdown(
            id='tissue-category-dropdown',
            options = [{'label' : k, 'value': k} for k in tcdict.keys()],
            placeholder = "Select a tissue category..."
        ),
    ],style={'width': '49%', 'display': 'inline-block', 'horizontalAlign': 'right'}),

    html.Div([
        dcc.Dropdown(
            id='tissue-dropdown',
            placeholder = "Select a tissue...",
            multi=True)

    ],style={'width': '49%', 'display': 'inline-block', 'horizontalAlign': 'right'}),
    
    html.Div([
        dcc.Graph(id='twas-graph', style={'width': '180vh', 'height': '45vh', 'display': 'inline-block'}),
        dcc.Graph(id='gwas-graph', style={'width': '180vh', 'height': '45vh', 'display': 'inline-block'}),
    ])
])

@app.callback(
    Output('tissue-dropdown', 'options'),
    [Input('tissue-category-dropdown', 'value')])
def set_tissue_options(selected_category):
    return [{'label': i, 'value': i} for i in tcdict[selected_category]]

@app.callback(
    Output('tissue-dropdown', 'value'),
    [Input('tissue-dropdown', 'options')])
def set_tissue_value(available_options):
    return available_options[0]['value']

@app.callback(Output('gwas-graph', 'figure'), [Input('interaction-dropdown', 'value'), Input('tissue-dropdown', 'value'), Input('maf-slider', 'value')]) #, [State('maf-slider', 'marks')])
def update_graph(selected_dropdown_value, selected_dropdown_value2, selected_maf):
    #Here it would be possible to split df into 9 separate traces and return the list
    df = d_order[(d_order.Interaction==selected_dropdown_value) & (d_order.Tissue=="GWAS") & (d_order.MAF > selected_maf[0]) & (d_order.MAF < selected_maf[1])]
    #cols = df['PHE'].map(colorsIdx)
    dl = []
    for i in df.PHE.unique():
        dfsub = df[df.PHE==i]
        #Text needs to be prefilled so add hover column here
        dff = {'x': dfsub.pos_index, 'y': dfsub.logp, 'mode': 'markers', 'text': dfsub.SNP, 'name': i}
        dl.append(dff)
    return {
        'data': dl,
        'layout': {'shapes': sl, 'title': 'GWAS Plot', 'showlegend': True, 'xaxis': {'range': [xmin, xmax], 'ticktext': lims.chrmin, 'tickvals': lims.av, 'showgrid': False}, 'yaxis': {'range': [ymin,ymax], 'title': '-log10(pvalue)'}, 'hovermode': 'closest'}
    }

@app.callback(Output('twas-graph', 'figure'), [Input('interaction-dropdown', 'value'), Input('tissue-dropdown', 'value'), Input('tissue-category-dropdown', 'value')])
def update_graph(selected_dropdown_value, selected_dropdown_value2, selected_dropdown_value3):
    #Here it would be possible to split df into 9 separate traces and return the list   
    df = d_order[(d_order.Interaction==selected_dropdown_value) & (d_order['Tissue'].isin(selected_dropdown_value2))]
    #df = d_order[(d_order.Interaction==selected_dropdown_value) & (d_order.Tissue==selected_dropdown_value2)]
    #cols = df['PHE'].map(colorsIdx)
    dl = []
    for i in df.PHE.unique():
        dfsub = df[df.PHE==i]
        #Text needs to be prefilled so add hover column here
        dff = {'x': dfsub.pos_index, 'y': dfsub.logp, 'mode': 'markers', 'text': dfsub.SNP, 'name': i}
        dl.append(dff)
    return {
        'data': dl,
        'layout': {'shapes': sl, 'title': 'TWAS Plot', 'showlegend': True, 'xaxis': {'range': [xmin, xmax], 'ticktext': lims.chrmin, 'tickvals': lims.av, 'showgrid': False}, 'yaxis': {'range': [ymin,ymax], 'title': '-log10(pvalue)'}, 'hovermode': 'closest'}
    }

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})

if __name__ == '__main__':
    app.run_server()
