#%%

import numpy as np
import pandas as pd
import plotly.express as px

#%%

df = pd.read_csv("output.dat",sep=";")
df_new = df[["x","Q2","z"]].copy()
# df_new[f'{df.columns[3]}/{df.columns[4]}']=df[df.columns[3]]/df[df.columns[4]]
df_new[f'{df.columns[5]}/{df.columns[6]}']=df[df.columns[5]]/df[df.columns[6]]
df_melted = df_new.melt(id_vars=["x", "Q2", "z"], var_name="function", value_name="value")

fig = px.line(
    df_melted,
    x="z",
    y="value",
    line_dash="function",
    animation_frame="x",
	color="Q2",
	log_x=True,
	markers=True
)

fig.update_layout(
    xaxis_title="z",
    legend_title="Function"
)

# sizes = [12,10,8,6]
# for i, trace in enumerate(fig.data):
# 	trace.marker.size = sizes[i]

fig.show()

#%%
import pandas as pd
import plotly.express as px
from dash import Dash, dcc, html, Input, Output

# --- Load and reshape your data ---
df = pd.read_csv("output.dat", sep=";")
df_melted = df_new.melt(id_vars=["x", "Q2", "z"], var_name="function", value_name="value")

# --- Create app ---
app = Dash(__name__)

app.layout = html.Div([
    html.H2("Interactive Plot with Two Sliders (x and Q2)"),

    # Sliders
    html.Div([
        html.Label("x value:"),
        dcc.Slider(
            id='x-slider',
            min=df_melted["x"].min(),
            max=df_melted["x"].max(),
            step=None,
            marks={float(v): str(v) for v in sorted(df_melted["x"].unique())},
            value=df_melted["x"].min()
        ),
        html.Label("Q2 value:"),
        dcc.Slider(
            id='q2-slider',
            min=df_melted["Q2"].min(),
            max=df_melted["Q2"].max(),
            step=None,
            marks={float(v): str(v) for v in sorted(df_melted["Q2"].unique())},
            value=df_melted["Q2"].min()
        )
    ], style={"width": "80%", "margin": "auto"}),

    # Plot
    dcc.Graph(id='line-plot')
])

# --- Define callback to update plot based on sliders ---
@app.callback(
    Output('line-plot', 'figure'),
    [Input('x-slider', 'value'),
     Input('q2-slider', 'value')]
)
def update_plot(selected_x, selected_q2):
    filtered = df_melted[(df_melted["x"] == selected_x) & (df_melted["Q2"] == selected_q2)]

    fig = px.line(
        filtered,
        x="z",
        y="value",
        color="function",
        markers=True,
        log_x=True,
        title=f"x = {selected_x}, Q² = {selected_q2}"
    )

    fig.update_layout(transition_duration=300)
    return fig

# --- Run app ---
if __name__ == "__main__":
    app.run(debug=True)
