#%%

import numpy as np
import pandas as pd
import plotly.express as px

#%%

df = pd.read_csv("output.dat",sep=";")
df_new = df[["x","Q2"]].copy()
cols=df.columns
df_new[f"{cols[2]}/{cols[3]}"]=df[cols[2]]/df[cols[3]]
df_new[f"{cols[4]}/{cols[5]}"]=df[cols[4]]/df[cols[5]]

df_melted = df_new.melt(id_vars=["x", "Q2"], var_name="function", value_name="value")
# df_melted = df.melt(id_vars=["x", "Q2"], var_name="function", value_name="value")
df_melted["chi"] = df_melted["Q2"]/(1.5**2)
df_melted["etamax"] = df_melted["Q2"]*(1./df_melted["x"]-1.)
df_melted["etamin"] = df_melted["Q2"]*(1./(1./(1.+4*1.5**2/df_melted["Q2"]))-1.)
df_melted["eta in range"] = (df_melted["etamax"] <= 1e6) * (df_melted["etamin"] >= 1e-6)
df_melted["chi in range"] = (df_melted["chi"] <= 1e5) * (df_melted["chi"] >= 1e-3)


# Plot with animation over y
fig = px.line(
    df_melted,
    x="x",
	# x="Q2",
    y="value",
    color="function",
    animation_frame="Q2",
	# animation_frame="x",
	hover_data=["chi","etamin","etamax", "chi in range","eta in range"],
	log_x=True,
	markers=True
)

fig.update_layout(
    xaxis_title="x",
	# xaxis_title="Q2",
    legend_title="Function",
	hoverlabel=dict(
    	bgcolor="white"  # 70% opacity
	)
)

sizes = [12,10,8,6]
for i, trace in enumerate(fig.data):
	trace.marker.size = sizes[i]

fig.show()
