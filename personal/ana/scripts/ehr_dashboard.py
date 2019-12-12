import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
from datetime import datetime as dt
import pandas as pd
import math
import numpy as np

icd = pd.read_csv("ICD_Cat.txt", sep="\t")
l = pd.read_csv("Labs.txt", sep="\t")
m = pd.read_csv("Meds_fmt.txt", sep="\t")

icd['SHIFTED_ENC_DT']= pd.to_datetime(icd['SHIFTED_ENC_DT'])
l['RESULT_DATE_SHIFT']= pd.to_datetime(l['RESULT_DATE_SHIFT'])
m['RX_DATE'] = pd.to_datetime(m['RX_DATE'])

#Set up the menus
app = dash.Dash(
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"}
    ]
)

ll = []
for i in l.DATASET.unique():
   choices = {'label': i, 'value': i}
   ll.append(choices)


ml = []
ml.append({'label': "All", 'value': "All"})
for i in m.Category.unique():
   choices = {'label': i, 'value': i}
   ml.append(choices)

app.layout = html.Div([

    html.Div([
        html.Label("Date Range"),
        dcc.DatePickerRange(
            id='date-range-picker',
            min_date_allowed=dt(2003, 1, 1),
            max_date_allowed=dt(2019, 4, 1),
            initial_visible_month=dt(2003, 1, 1),
            #start_date=dt(2003, 1, 1),
            end_date=dt(2019, 4, 1),
        )
    ],style={'height': '60px', 'width': '33%', 'display': 'inline-block', 'verticalAlign': 'top', 'horizontalAlign': 'right', 'backgroundColor':'#E8E8E8'}),

    html.Div([
        html.Label('Labs'),
        dcc.Dropdown(
            id='lab-dropdown',
            options = ll,
            placeholder = "Select a lab...",
            multi=True
        ),
    ],style={'width': '33%', 'display': 'inline-block', 'verticalAlign': 'top', 'backgroundColor':'#E8E8E8'}),

    html.Div([
        html.Label('Medications'),
        dcc.Dropdown(
            id='med-dropdown',
            options = ml,
            placeholder = "Select a medicationm category...",
            multi=True
        ),
    ],style={'width': '33%', 'display': 'inline-block', 'verticalAlign': 'top', 'horizontalAlign': 'left', 'backgroundColor':'#E8E8E8'}),

    html.Div([
        dcc.Graph(id='icd-graph', style={'width': '98vw', 'height': '28vh', 'display': 'inline-block', 'backgroundColor':'#E8E8E8'}),
#    ]),

  #  html.Div([
        dcc.Graph(id='lab-graph', style={'width': '98vw', 'height': '28vh', 'display': 'inline-block', 'backgroundColor':'#E8E8E8'}),
 #   ]),

 #   html.Div([
        dcc.Graph(id='med-graph', style={'width': '98vw', 'height': '28vh', 'display': 'inline-block', 'backgroundColor':'#E8E8E8'})
    ])

], style={'backgroundColor':'#E8E8E8'})

@app.callback(Output('icd-graph', 'figure'), [Input('date-range-picker', 'start_date'), Input('date-range-picker', 'end_date')])
def update_graph(start_date, end_date):
    dfl = icd[(icd.SHIFTED_ENC_DT >= start_date) & (icd.SHIFTED_ENC_DT <= end_date)]
    dcl = []
    for i in icd.Category.unique():
        dfsub = dfl[dfl.Category==i]
        dfsub['Hover'] = dfsub.CODE
        #Text needs to be prefilled so add hover column here
        dff = {'x': dfsub.SHIFTED_ENC_DT, 'y': dfsub.PATIENT_CLASS, 'mode': 'markers', 'text': dfsub.Hover, 'name': i}
        dcl.append(dff)
    return {
        'data': dcl,
        'layout': {'showlegend': True, 'xaxis': {'title':'ICD Codes','range': [start_date, end_date]}, 'hovermode': 'closest', 'margin':{'t':15}},
        'style': {'marginBottom': '0.1em', 'marginTop': '0.1em'}
    }

@app.callback(Output('lab-graph', 'figure'), [Input('date-range-picker', 'start_date'), Input('date-range-picker', 'end_date'), Input('lab-dropdown', 'value')])
def update_graph(start_date, end_date, value):
    dfl = l[(l.RESULT_DATE_SHIFT >= start_date) & (l.RESULT_DATE_SHIFT <= end_date) & (l['DATASET'].isin(value))]
    dcl = []
    for i in dfl.DATASET.unique():
        dfsub = dfl[dfl.DATASET==i]
        dfsub['Hover'] = dfsub.ABNORMAL
        dfs = dfsub.sort_values(by='RESULT_DATE_SHIFT', ascending=True)
        dff = {'x': dfs.RESULT_DATE_SHIFT, 'y': dfs.RESULT_VALUE_NUM, 'mode': 'lines+markers', 'hover': dfs.Hover, 'name': i}
        dcl.append(dff)
    return {
        'data': dcl,
        'layout': {'showlegend': True, 'xaxis': {'title': 'Clinical Labs', 'range': [start_date, end_date]}, 'hovermode': 'closest', 'margin':{'t':0}}
    }

@app.callback(Output('med-graph', 'figure'), [Input('date-range-picker', 'start_date'), Input('date-range-picker', 'end_date'), Input('med-dropdown', 'value')])
def update_graph(start_date, end_date, value):
    dfl = m[(m.RX_DATE >= start_date) & (m.RX_DATE <= end_date) & (m['Category'].isin(value))]
    dcl = []
    for i in dfl.Medication.unique():
        dfsub = dfl[dfl.Medication==i]
        dfsub['Hover'] = dfsub.Category
        dff = {'x': dfsub.RX_DATE, 'y': dfsub.Medication, 'mode': 'markers', 'hover': dfsub.Hover, 'name': i}
        dcl.append(dff)
    return {
        'data': dcl,
        'layout': {'showlegend': True, 'xaxis': {'title': 'Medications','range': [start_date, end_date]}, 'hovermode': 'closest', 'margin':{'t':0}}
    }


app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})

if __name__ == '__main__':
    app.run_server()
