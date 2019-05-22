import base64
import io

import plotly.graph_objs as go

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table

import pandas as pd
import numpy as np

from scipy import optimize

#Define lognormal function
def lognormal(x, mu, s):
	return (1/(x*s*np.sqrt(2*np.pi))) * (np.exp(-(((np.log(x/mu))**2)/(2*s**2))))


#Apply style rules which are provided by plotly
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

#Style variables
button_font='22px'
title_graph_style=dict(family='Arial, sans-serif', size=22)
axis_graph_style=dict(family='Arial, sans-serif', size=18, color='#7f7f7f')

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
#Create the data structures which will contain the data
AbsorptionSpectrum = pd.DataFrame()
AbsorptionDatabase = pd.DataFrame()
Jacobian = pd.DataFrame()

empty=[]

NPsizes = np.array(empty) #create an array for the NPs sizes
wavelength = np.array(empty) #create an array for the wavelength
A = np.array(empty) #create an array for the Absorption Database
b = np.array(empty) #create an array for the Absorption Spectra
NPsizes_frecuency = np.array(empty) #create an array for the frecuency of Sizes

#Keeps track of the state of the system and defines the order in which the program is operated
flag = 0

#Activate the filter
flag_filter = 0

#Threshold filter
aux_threshold = 0

#Set layout
app.layout = html.Div(
	html.Div([
		#Title
		html.Div(
			[
				html.Img(
					src="http://www.dbbe.fcen.uba.ar/contenido/objetos/mini_100/LogoExactasNEGRO.jpg",
					className='three columns',
					 style={
						'height': '8%',
						'width': '8%',
						'float': 'left',
						'position': 'center',
						'margin-top': 5,
					},
				),
				html.Div([
				html.H1('Diameter Distribution by Deconvolution'),
				html.H4('Web app for determining the size of nanoparticles from absorption spectra. DOI:10.XXXX/XXXXXXX'),
				],
				className='nine columns',
				style={
				'textAlign':'center',
				},
				),
				html.Img(
					src="http://www.dbbe.fcen.uba.ar/contenido/objetos/mini_100/LogoExactasNEGRO.jpg",
					className='three columns',
					style={
						'height': '8%',
						'width': '8%',
						'float': 'right',
						'position': 'center',
						'margin-top': 5,
					},
				),
			], className="row",
		),

		html.Hr(),  # horizontal line


		#Load Jacobian
		 html.Div(
			[
				html.Div(
					[
						dcc.Upload(
						id='upload-J',
						children=html.Button('(4) Jacobian', id='buttonJ', style={'font-size' : button_font}),
						# Allow multiple files to be uploaded
						multiple=True
						),
						dcc.Checklist(
							id='filterChecklist',
							options=[
								{'label': 'Filter the results', 'value': 1},
							],
							values=[],
							style={'font-size' : '22px'}
						),
						dcc.Input(
							id='threshold',
							placeholder='Threshold filter...',
							step=0.1,
							type='number',
							value='',
							style={'font-size' : '22px'}
						),
						html.Button('Apply filter', id='buttonFilter', style={'font-size' : '18px'}),
					], style={'textAlign': 'center'},
				),
			], className="row"
		),



		#PSD Graph
		html.Div(
			id='graph-PSD',
			className="row"
		),


		html.Hr(),  # horizontal line


		#Calculate DDD
		 html.Div(
			[
				html.Div(
					[
						html.Button('(3) Execute NNLS', id='buttonNNLS', style={'font-size' : button_font}),
					], style={'textAlign': 'center'},
				)
			], className="row",
		),


		#DDD Graph
		html.Div([
			html.Div(
				id='graph-DDD',
			),
		],className="row"
		),


		html.Hr(),  # horizontal line

		#Load Absorption Data
		 html.Div(
			[
				html.Div(
					[
						dcc.Upload(
						id='upload-AS',
						children=html.Button('(1) Absorption Spectrum', id='buttonAS', style={'font-size' : button_font}),
						# Allow multiple files to be uploaded
						multiple=True
						),
					], className= 'six columns', style={'textAlign': 'center'},
				),

				html.Div(
					[
						dcc.Upload(
						id='upload-AD',
						children=html.Button('(2) Absorption Database', id='buttonAD', style={'font-size' : button_font}),
						# Allow multiple files to be uploaded
						multiple=True
						),
					], className= 'six columns', style={'textAlign': 'center'},
				)
			], className="row"
		),

		#Absorption Spectrum & Absorption Database Graph
		html.Div(
			[
				html.Div(
					id='graph-AS', className= 'six columns'
				),

				html.Div(
					id='graph-AD', className= 'six columns'
					)
			], className="row"
		),



		html.Hr(),  # horizontal line

		html.Div(
			id='hidden',
			className="row"
		),

		#Author credits
		html.Div([
		dcc.Markdown(''' ### Diego Onna [diego.onna@qi.fcen.uba.ar](mailto:diego.onna@qi.fcen.uba.ar)''')
		],
		className="row",
		style={'textAlign': 'center'}),

		# Read me
		html.Div([
			dcc.Markdown('''### Instructions''')
		],
			className="row",
			style={'textAlign': 'center'}),
		html.Div([
			dcc.Markdown('''##### 1) Load the Absorption Spectrum\n ##### 2) Load the Absorption Database\n ##### 3) Execute NNLS\n ##### 4) Load the Jacobian\n ##### 5) Apply a threshold filter for eliminating the numeric error in lower sizes			
			''')
		],
			className="row"),

	], className='ten columns offset-by-one')
)

#Load Data & Plot Absorption Spectrum
def plotAS(contents, filename):
	global AbsorptionSpectrum
	global flag
	flag = 1
	content_type, content_string = contents.split(',')
	decoded = base64.b64decode(content_string)
	try:
		if 'csv' in filename:
			# Assume that the user uploaded a CSV file
			AbsorptionSpectrum = pd.read_csv(
				io.StringIO(decoded.decode('utf-8')))
		elif 'xls' in filename:
			# Assume that the user uploaded an excel file
			AbsorptionSpectrum = pd.read_excel(io.BytesIO(decoded))
	except Exception as e:
		print(e)
		return html.Div([
			'There was an error processing this file.'
		])

	return html.Div([
				dcc.Graph(
					figure={
						'data': [
						   {'x': AbsorptionSpectrum.iloc[:,0], 'y': AbsorptionSpectrum.iloc[:,1], 'type': 'line'},
						],
						'layout': {
							'title': 'Absorption Spectrum',
                            'titlefont' : title_graph_style,
						'xaxis': {
							'title': 'Wavelength (nm)',
                            'titlefont' : axis_graph_style
						},
						'yaxis': {
							'title': 'Absorption',
                            'titlefont' : axis_graph_style
						},
						}
					}
				)
			]
	)


#Callback Upload Absorption Spectra
@app.callback(Output('graph-AS', 'children'),
			 [Input('upload-AS', 'contents')],
			 [State('upload-AS', 'filename')])
def update_graphES(list_of_contents, list_of_names):
	if list_of_contents is not None:
		children = [
			plotAS(c, n) for c, n in
			zip(list_of_contents, list_of_names)]
		return children



#Load Data & Plot Absorption Database
def plotAD(contents, filename):
	global AbsorptionDatabase
	global flag

	if (flag == 1):
		flag = 2

		content_type, content_string = contents.split(',')
		decoded = base64.b64decode(content_string)
		try:
			if 'csv' in filename:
				# Assume that the user uploaded a CSV file
				AbsorptionDatabase = pd.read_csv(
					io.StringIO(decoded.decode('utf-8')))
			elif 'xls' in filename:
				# Assume that the user uploaded an excel file
				AbsorptionDatabase = pd.read_excel(io.BytesIO(decoded))
		except Exception as e:
			print(e)
			return html.Div([
				'There was an error processing this file.'
			])


		columns = AbsorptionDatabase.columns.values.tolist()
		columns = columns[1:]

		trace=[go.Scatter(x=AbsorptionDatabase.iloc[:,0], y = AbsorptionDatabase[i], mode = 'lines', name = i) for i in columns]


		return html.Div([
					dcc.Graph(
						figure={
							'data': trace,
							'layout': {
								'title': 'Absorption Database',
                                'titlefont': title_graph_style,
                                'xaxis': {
									'title': 'Wavelength (nm)',
                                    'titlefont': axis_graph_style
                                },
								'yaxis': {
									'title': 'Absorption',
                                    'titlefont': axis_graph_style
                                },
							}
						}
					)
				]
		)
	else:
		return html.Div([
				html.H1('Load first the Absorption Spectrum.',
				style={'textAlign': 'center'},
				),
			])


#Callback Upload Absorption Database
@app.callback(Output('graph-AD', 'children'),
			 [Input('upload-AD', 'contents')],
			 [State('upload-AD', 'filename')])
def update_graphED(list_of_contents, list_of_names):
	if list_of_contents is not None:
		children = [
			plotAD(c, n) for c, n in
			zip(list_of_contents, list_of_names)]
		return children


#Callback make DdD calculation
@app.callback(Output('graph-DDD', 'children'),
			 [Input('buttonNNLS', 'n_clicks')])
def calculate_DdD(n_clicks):
	global flag
	global NPsizes
	global wavelength
	global A
	global b
	global AbsorptionDatabase
	global NPsizes_frecuency

	if n_clicks != None:
		if (flag == 2):
			flag = 3

			NPsizes = np.array([float(i) for i in AbsorptionDatabase.columns.values[1:]]) #define the array with the NPs sizes
			wavelength=AbsorptionSpectrum.values[:,0].flatten() #define the array with the wavelength

			A=AbsorptionDatabase.values[:,1:] #define the array with the Absorption Database
			b=AbsorptionSpectrum.values[:,1:].flatten() #define the array with the Absorption Spectra and make it 1D

			#Find x which solves Ax=b where x is NPsizes_frecuency

			NPsizes_frecuency, rnorm = optimize.nnls(A,b)

			dataToFit=go.Scatter( x = wavelength,
								  y = b,
								  mode = 'lines',
								  name = 'Data'
								)

			dataFitted=go.Scatter( x = wavelength,
								   y = np.matmul(A,NPsizes_frecuency),
								   mode = 'lines',
								   name = 'Fit'
								)

			columns = AbsorptionDatabase.columns.values.tolist()
			columns = columns[1:]

			trace=[go.Scatter(x=AbsorptionDatabase.iloc[:,0], y = AbsorptionDatabase[i]*NPsizes_frecuency[columns.index(i)], mode = 'lines', name = i) for i in columns]
			trace.append(dataToFit)
			trace.append(dataFitted)


			return html.Div(
				[
					dcc.Graph(
						figure={
							'data': trace,
							'layout': {
								'title': 'Absorption Spectrum Fitted',
								'titlefont': title_graph_style,
								'xaxis': {
									'title': 'Wavelength (nm)',
									'titlefont': axis_graph_style
								},
								'yaxis': {
									'title': 'Absorption',
									'titlefont': axis_graph_style
								},
								}
						}
					)
				],className='eight columns offset-by-two',
			)
		else:
			return html.Div([
				html.H1('Load first the Absorption Spectrum, and then the Absorption Database.',
				style={'textAlign': 'center'},
				),
			])


#Load Jacobian & Plot PSD
def plotPSD(contents, filename):
	global Jacobian
	global NPsizes_frecuency
	global flag

	if (flag == 3):
		flag = 4

		content_type, content_string = contents.split(',')
		decoded = base64.b64decode(content_string)
		try:
			if 'csv' in filename:
				# Assume that the user uploaded a CSV file
				Jacobian = pd.read_csv(
					io.StringIO(decoded.decode('utf-8')))
			elif 'xls' in filename:
				# Assume that the user uploaded an excel file
				Jacobian = pd.read_excel(io.BytesIO(decoded))
		except Exception as e:
			print(e)
			return html.Div([

				'There was an error processing this file.'
			])

		FitXvalues = np.linspace(min(Jacobian['Size']), max(Jacobian['Size']), 50)
		x_data=Jacobian['Size']
		y_data=NPsizes_frecuency*Jacobian['J']/max(NPsizes_frecuency*Jacobian['J'])

		params, params_covariance = optimize.curve_fit(lognormal, x_data, y_data)

		return html.Div([
					dcc.Graph(
						figure={
							'data': [
										go.Bar(
												x=x_data,
												y=y_data,
												name='DDD'
										),
										go.Scatter(
											x = FitXvalues,
											y = lognormal(FitXvalues, params[0], params[1]),
											name='lognormal fitting'
										),
							],
							'layout': {
								'title': 'Particle Size Distribution',
								'titlefont': title_graph_style,
								'xaxis': {
									'title': 'Particle Size (nm)',
									'titlefont': axis_graph_style
								},
								'yaxis': {
									'title': 'Frequency',
									'titlefont': axis_graph_style
								},
								'annotations':
										[
										{
											'x': 5,
											'y': 0.9,
											'align': "center",
											'arrowcolor': "rgba(0, 0, 0, 0)",
											'font': {'size': 18},
											'text': 'Mean Size = {0:.1f} nm   Desviation = {1:.1f} nm'.format(np.exp(np.log(params[0])+0.5*params[1]*params[1]),np.exp(np.log(params[0])+0.5*params[1]*params[1])*np.sqrt(np.exp(params[1]*params[1])-1)),
											'textangle': 0
										}
									]
							}
						}
					)
				],className='eight columns offset-by-two',
		)
	else:
		return html.Div([
				html.H1('Load first the Absorption Spectrum, and then the Absorption Database. Then calculate the DdD',
				style={'textAlign': 'center'},
				),

			])

#Apply filter
def applyFilter(n_clicks):
	global Jacobian
	global NPsizes_frecuency
	global flag
	global flag_filter
	global aux_threshold

	if (flag == 4):
		#flag = 5

		FitXvalues = np.linspace(min(Jacobian['Size']), max(Jacobian['Size']), 50)
		x_data=Jacobian['Size']
		y_data=NPsizes_frecuency*Jacobian['J']/max(NPsizes_frecuency*Jacobian['J'])


		if (flag_filter == 1):
			for i in x_data.index:
				if x_data.iloc[i] < aux_threshold:
					y_data.iloc[i] = 0

		params, params_covariance = optimize.curve_fit(lognormal, x_data, y_data)

		return html.Div([
					dcc.Graph(
						figure={
							'data': [
										go.Bar(
												x=x_data,
												y=y_data,
												name='DDD'
										),
										go.Scatter(
											x = FitXvalues,
											y = lognormal(FitXvalues, params[0], params[1]),
											name='lognormal fitting'
										),
							],
							'layout': {
								'title': 'Particle Size Distribution',
								'titlefont': title_graph_style,
								'xaxis': {
									'title': 'Particle Size (nm)',
									'titlefont': axis_graph_style
								},
								'yaxis': {
									'title': 'Frequency',
									'titlefont': axis_graph_style
								},
								'annotations':
										[
										{
											'x': 5, 
											'y': 0.9, 
											'align': "center",
											'arrowcolor': "rgba(0, 0, 0, 0)", 
											'font': {'size': 18},
											'text': 'Mean Size = {0:.1f} nm   Desviation = {1:.1f} nm'.format(np.exp(np.log(params[0])+0.5*params[1]*params[1]),np.exp(np.log(params[0])+0.5*params[1]*params[1])*np.sqrt(np.exp(params[1]*params[1])-1)), 
											'textangle': 0
										}
									]                                
							}
						}
					)
				],className='eight columns offset-by-two',
		)
	else:
		return html.Div([
				html.H1('Load first the Absorption Spectrum, and then the Absorption Database. Then calculate the DdD. Finally load the Jacobian',
				style={'textAlign': 'center'},
				),

			])


#Callback Upload Jacobian
@app.callback(Output('graph-PSD', 'children'),
			 [Input('upload-J', 'contents'),Input('buttonFilter', 'n_clicks')],
			 [State('upload-J', 'filename')])
def update_graphPSD(list_of_contents, n_clicks, list_of_names): #The order of the arguments is defined by the order of the callbacks
	children=None
	if (list_of_contents is not None):
		children = [
			plotPSD(c, n) for c, n in
			zip(list_of_contents, list_of_names)]
	if n_clicks != None:
		children = applyFilter(n_clicks)
	return children

#Callback Upload Filter
@app.callback(Output('hidden', 'children'),
			 [Input('filterChecklist', 'values'),
			  Input('threshold','value')])
def update_filter(filterChecklist_values,threshold_value):
	global flag_filter
	global aux_threshold
	if len(filterChecklist_values) == 1:
		flag_filter = 1
		aux_threshold = threshold_value
	else:
		flag_filter = 0
		aux_threshold = 0
	pass



if __name__ == "__main__":
	app.run_server(debug=True)
