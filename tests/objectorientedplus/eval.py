#%%

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "vscode"
import io
import numpy as np

#%%

df = pd.read_csv("testdata.dat",sep=";")
fig = px.scatter(
	df,
	x="err",
	y="runtime",
	color="algorithm",
)

fig.update_layout(
    xaxis=dict(
        type="log",
        # tickformat=".0e",
		exponentformat="power",
        title="integration error"
    ),
    yaxis=dict(
        type="log",
        # tickformat=".0e",
		exponentformat="power",
        title="runtime(ns)"
    )
)

fig.show()

#%%

exact_result = 2

df1 = pd.read_csv("integrate_gauss1.dat",sep=";",skiprows=1,names=["runtime(ns)","result"])
df2 = pd.read_csv("integrate_gsl.dat",sep=";",skiprows=1,names=["runtime(ns)","result"])

fig = go.Figure()
fig.add_trace(go.Scatter(y=df1["runtime(ns)"], x=np.abs(exact_result-df1["result"]), mode="markers", name="gauss1 (openQCDrad)"))
fig.add_trace(go.Scatter(y=df2["runtime(ns)"], x=np.abs(exact_result-df2["result"]), mode="markers", name="gsl_integrate"))

fig.update_layout(
    xaxis=dict(
        type="log",
        # tickformat=".0e",
		exponentformat="power",
        title="integration error"
    ),
    yaxis=dict(
        type="log",
        # tickformat=".0e",
		exponentformat="power",
        title="runtime(ns)"
    ),
	# title=r'$\int_0^1dx \exp(Ax)*\cos(Bx)$'
	title=r'$\int_0^1dx \sqrt{x}$'
)


fig.show()
