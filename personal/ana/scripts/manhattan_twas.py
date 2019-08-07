import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import math
import numpy as np

###Formatting
dat = pd.read_csv("~/Desktop/projects/plotly_dash/MetaXcan_GERA_GLGC_eMERGEIII_UKBB_v7_0.01.txt", sep="\t")
dat['chromosome'] = dat.chromosome.astype(int)
dat['gene_start_position'] = dat.gene_start_position.astype(int)
dat['posterior_inclusion'] = dat.posterior_inclusion.astype(float)
d_order = dat.sort_values(by=['chromosome','gene_start_position'])
d_order['pos_index'] = range(1, len(d_order)+1)

#Get tick positions
maxRows = d_order.loc[d_order.groupby('chromosome').pos_index.agg('idxmax')][['chromosome', 'gene_start_position', 'pos_index']]
minRows = d_order.loc[d_order.groupby('chromosome').pos_index.agg('idxmin')][['chromosome', 'gene_start_position', 'pos_index']]
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
d_order['Hover'] = "Gene: " + d_order['ENSG_gene'] + "\npvalue: " + d_order['pvalue'].astype(str) + "\nDataset: " + d_order['data'] + "\nTissue: " + d_order['tissue']

#Create shapes for chromosomes
sl = []
for i in range(int(len(limsshape.chrmin.unique()))):
    shapes = {'type': 'rect', 'x0': limsshape.iloc[i]['posindmin'], 'x1': limsshape.iloc[i]['posindmax'], 'y0': ymin-0.5, 'y1': ymax, 'fillcolor': '#d3d3d3', 'opacity': 0.3, 'layer': 'below', 'line':{'width': 0}}
    sl.append(shapes)

#Set up the menus
app = dash.Dash(
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"}
    ]
)
#Make these into a loop to get the label and values

#Get list of dataset terms for menu
dl = []
for i in d_order.data.unique():
   datasets = {'label': i, 'value': i}
   dl.append(datasets)

#Get list of tissue terms for menu
tl = []
for i in d_order.tissue.unique():
        tissues = {'label': i, 'value': i}
        tl.append(tissues)


#Get marker positions for slider based on max posterior inclusion
markers={0: str(0), 0.1: str(0.1), 0.2: str(0.2), 0.3: str(0.3), 0.4: str(0.4), 0.5: str(0.5), 0.6: str(0.6), 0.7: str(0.7), 0.8: str(0.8), 0.9: str(0.9), 1: str(1)}
maxpostinc = 1

app.layout = html.Div([
    html.Div([
        dcc.Dropdown(
            id='dataset-dropdown',
            options = dl,
            placeholder = "Select an dataset..."
        ),
    ],style={'width': '49%', 'display': 'inline-block'}),

    html.Div([
        dcc.RangeSlider(
            id='postinc-slider',
            min=0,
            max=maxpostinc,
            step=maxpostinc/10,
            marks=markers,
            value=[0,maxpostinc],
            allowCross=False
        )
    ],style={'width': '49%', 'display': 'inline-block', 'verticalAlign': 'top', 'horizontalAlign': 'right'}),

    html.Div([
        dcc.Dropdown(
            id='tissue-dropdown',
            options = tl,
            placeholder = "Select a tissue..."
        ),
    ],style={'width': '49%', 'display': 'inline-block', 'horizontalAlign': 'right'}),

    html.Div([
        dcc.Graph(id='twas-graph', style={'width': '180vh', 'height': '45vh', 'display': 'inline-block'}),
    ])
])

@app.callback(Output('twas-graph', 'figure'), [Input('dataset-dropdown', 'value'), Input('tissue-dropdown', 'value'), Input('postinc-slider', 'value')])
def update_graph(selected_dropdown_value, selected_dropdown_value2, selected_postinc):
    #Here it would be possible to split df into 9 separate traces and return the list   
    df = d_order[(d_order.data==selected_dropdown_value) & (d_order.tissue==selected_dropdown_value2) & (d_order.posterior_inclusion > selected_postinc[0]) & (d_order.posterior_inclusion < selected_postinc[1])]
    dsl = []
    for i in df.trait.unique():
        dfsub = df[df.trait==i]
        #Text needs to be prefilled so add hover column here
        dff = {'x': dfsub.pos_index, 'y': dfsub.logp, 'mode': 'markers', 'text': dfsub.Hover, 'name': i}
        dsl.append(dff)
    return {
        'data': dsl,
        'layout': {'shapes': sl, 'title': 'TWAS Plot', 'showlegend': True, 'xaxis': {'range': [xmin, xmax], 'ticktext': lims.chrmin, 'tickvals': lims.av, 'showgrid': False}, 'yaxis': {'range': [ymin,ymax], 'title': '-log10(pvalue)'}, 'hovermode': 'closest'}
    }

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})

if __name__ == '__main__':
    app.run_server()
