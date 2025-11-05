#%%

import numpy as np
import pandas as pd
import plotly.express as px
import glob

#%%
files = glob.glob("f2*.dat")
if(files.__len__ == 0):
	quit()
files.sort()
df = pd.read_csv(files[-1],sep=";",comment="#")

for col0 in [2,5,8]:
	df_new = df[["x","Q2"]].copy()
	for col in df.columns[col0:col0+3]:
		df_new[f'{col}/{df.columns[col0+2]}'] = df[col]/df[df.columns[col0+2]]
	df_melted = df_new.melt(id_vars=["x", "Q2"], var_name="function", value_name="value")
	sizes = [8,6,4]
	sizemap = dict(zip(df_melted["function"].unique(), sizes))
	df_melted["marker_size"] = df_melted["function"].map(sizemap)

	fig = px.line(
		df_melted,
		x="x",
		# x="Q2",
		y="value",
		color="function",
		animation_frame="Q2",
		# animation_frame="x",
		log_x=True,
		markers=True
	)

	fig.update_layout(
		xaxis_title="x",
		# xaxis_title="Q2",
		legend_title="Function"
	)

	for i, trace in enumerate(fig.data):
		trace.marker.size = sizes[i]

	fig.show()

#%%
files = glob.glob("pdf*.dat")
if(files.__len__ == 0):
	quit()
files.sort()
df = pd.read_csv(files[-1],sep=";",comment="#")

for col0 in [2,5,8,11,14,17,20,23,26,29,32,35,38]:
	df_new = df[["x","Q2"]].copy()
	for col in df.columns[col0:col0+3]:
		df_new[f'({col}-{df.columns[col0+2]})/{df.columns[col0+2]}'] = (df[col]-df[df.columns[col0+2]])/df[df.columns[col0+2]]
	df_melted = df_new.melt(id_vars=["x", "Q2"], var_name="function", value_name="value")
	sizes = [8,6,4]
	sizemap = dict(zip(df_melted["function"].unique(), sizes))
	df_melted["marker_size"] = df_melted["function"].map(sizemap)

	fig = px.line(
		df_melted,
		x="x",
		# x="Q2",
		y="value",
		color="function",
		animation_frame="Q2",
		# animation_frame="x",
		log_x=True,
		markers=True
	)

	fig.update_layout(
		xaxis_title="x",
		# xaxis_title="Q2",
		legend_title="Function"
	)

	for i, trace in enumerate(fig.data):
		trace.marker.size = sizes[i]

	fig.show()

#%%
files = glob.glob("alphas*.dat")
if(files.__len__ == 0):
	quit()
files.sort()
df = pd.read_csv(files[-1],sep=";",comment="#")

df_new = df[["Q2"]].copy()
for col in df.columns[1:]:
	df_new[f'{col}-{df.columns[1]}'] = df[col]-df[df.columns[1]]
df_melted = df_new.melt(id_vars=["Q2"],var_name="function",value_name="value")
sizes = [8,6,4]
sizemap = dict(zip(df_melted["function"].unique(), sizes))

fig = px.line(
	df_melted,
	x="Q2",
	y="value",
	color="function",
	log_x=True,
	markers=True
)

for i, trace in enumerate(fig.data):
	trace.marker.size = sizes[i]

fig.show()