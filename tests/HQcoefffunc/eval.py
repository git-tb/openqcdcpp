#%%

import numpy as np
import pandas as pd
import plotly.express as px

#%%

df = pd.read_csv("output.dat",sep=";",header=None)
i0 = 7

df_new = df[[0,1]].copy()

# df_new[i0] = df[i0]
# df_new[i0+10] = df[i0+10]
df_new[f'{i0}/{(i0+10)}'] = df[i0]/df[i0+10]

# df_new = df_new[df_new[0] < 2e-5]
df_melted = df_new.melt(id_vars=[0,1], var_name="function", value_name="value")


# Plot with animation over y
fig = px.line(
    df_melted,
    x=0,
    y="value",
    color="function",
    animation_frame=1,
	log_x=True,
	markers=True
)

fig.update_layout(
    xaxis_title="eta",
    legend_title="Function"
)

Ndata = len(fig.data)
sizes = 5+3*(Ndata-np.arange(Ndata))
for i, trace in enumerate(fig.data):
	trace.marker.size = sizes[i]

fig.show()
