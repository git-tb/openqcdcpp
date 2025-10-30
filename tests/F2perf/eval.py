#%%

import numpy as np
import pandas as pd
import plotly.express as px

#%%

df = pd.read_csv("output.dat",sep=";")
df_melted = df.melt(id_vars=["x", "Q2"], var_name="function", value_name="value")

# Plot with animation over y
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

fig.show()
# %%

df_nnlo = df[["x","Q2"]].copy()
df_nnlo = df_nnlo.assign(newcol1 = df["F2@nnlo(cpp,apr1)"]/df["F2@nnlo(cpp,ex)"])
df_nnlo = df_nnlo.assign(newcol2 = df["F2@nnlo(cpp,apr2)"]/df["F2@nnlo(cpp,ex)"])
df_nnlo = df_nnlo.assign(newcol3 = df["F2@nnlo(cpp,ex)"]/df["F2@nnlo(cpp,ex)"])

df_nnlo = df_nnlo.rename(columns={
		"newcol1":"F2@nnlo(apr1)/F2@nnlo(ex)",
		"newcol2":"F2@nnlo(apr2)/F2@nnlo(ex)",
		"newcol3":"F2@nnlo(ex)/F2@nnlo(ex)",
		})

df_nnlo_melted = df_nnlo.melt(id_vars=["x", "Q2"], var_name="function", value_name="value")

# Plot with animation over y
fig = px.line(
    df_nnlo_melted,
    x="x",
    y="value",
    color="function",
    animation_frame="Q2",
	log_x=True,
	markers=True
)

fig.update_layout(
    xaxis_title="x",
    legend_title="Function"
)

fig.show()
