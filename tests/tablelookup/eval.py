#%%
import pandas as pd
import plotly.express as px

#%%

df = pd.read_csv("test.dat",sep=";")
px.scatter(
	df,
	x="listlen",
	y=["time algo1","time algo2"],
	log_x=True,
	log_y=True
)