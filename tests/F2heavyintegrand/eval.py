#%%

import numpy as np
import pandas as pd
import plotly.express as px

#%%

df = pd.read_csv("output.dat",sep=";")
df_melted = df.melt(id_vars=["x", "Q2", "z"], var_name="function", value_name="value")


# Plot with animation over y
fig = px.line(
    df_melted,
    x="z",
    y="value",
    color="function",
    animation_frame="x",
	log_x=True,
	markers=True
)

fig.update_layout(
    xaxis_title="z",
    legend_title="Function"
)

sizes = [12,10,8,6]
for i, trace in enumerate(fig.data):
	trace.marker.size = sizes[i]

fig.show()
