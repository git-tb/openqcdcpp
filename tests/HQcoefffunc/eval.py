#%%

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

#%%
### NOTE: for the gluon case, use this evaluation cell. For the quark case, use the next



df = pd.read_csv("output.dat",sep=";",header=None)
i0 = 12

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


#%%

df = pd.read_csv("output.dat",sep=";",header=None)
i0 = 2

df_new = df[[0,1]].copy()

# df_new[i0] = df[i0]
# df_new[i0+12] = df[i0+12]
df_new[f'f{i0}/f{(i0+12)}'] = df[i0]/df[i0+12]
# df_new[f'{i0}'] = df[i0]

### finding the zeros crossings in x0 (==eta) for fixed x1 == (chi), I think this
### is the way to do it
arr = np.reshape(df[i0+12], (np.unique(df[1]).size, np.unique(df[0]).size))
arr=np.vstack([arr,arr[-1]])
zeros = 1/(np.abs(np.diff(np.sign(arr),axis=0))/2)
df_new[f'zeros of f{i0+12}'] = zeros.flatten()

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