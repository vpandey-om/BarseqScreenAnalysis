import pandas as pd
import numpy as np
from scipy.stats import lognorm
import plotly.graph_objects as go
import pickle

# mu, sigma = 3., 1. # mean and standard deviation
# x01=np.array([])
# for i in range(4):
#     x0 = np.random.lognormal(mu, sigma, 63)
#     x0=np.log2(x0/x0.sum())
#     x01=np.concatenate((x01, x0), axis=None)

pool1=pickle.load(open('/Users/vikash/git-hub/BarseqScreenAnalysis/rel_freq_pool1.pickle','rb'))
pilot1=pickle.load(open('/Users/vikash/git-hub/BarseqScreenAnalysis/rel_freq_pilot.pickle','rb'))

# x1 = np.random.lognormal(mu, sigma, 120)
# x1=np.log2(x1/x1.sum())
# x21=np.array([])
# for i in range(12):
#     x2 = np.random.lognormal(mu, sigma, 202)
#     x2=np.log2(x2/x2.sum())
#     x21=np.concatenate((x21, x2), axis=None)

x01= pilot1['d0_mean'].to_numpy()
x21=pool1['d0_mean'].to_numpy()


print (x01.mean(),x21.mean())


fig = go.Figure()
# fig.add_trace(go.Histogram(x=x01,name = "simulated pilot distribution (%d)"%len(x01)))
fig.add_trace(go.Histogram(x=x01,name = "pilot distribution (%d)"%len(x01)))
trace_x0 = go.Scatter(
       x = [x01.mean(),x01.mean()],
       y = [0,40],
       mode = "lines",
       name = "pilot_mean",
       )

fig.add_trace(trace_x0)
# fig.add_trace(go.Histogram(x=x1))
# trace_x1 = go.Scatter(
#        x = [x1.mean(),x1.mean()],
#        y = [0,40],
#        mode = "lines",
#        name = "x1_mean",
#        )
#
# fig.add_trace(trace_x1)
# fig.add_trace(go.Histogram(x=x21,name = "simulated pool1 distribution (%d)"%len(x21)))
fig.add_trace(go.Histogram(x=x21,name = "pool1 distribution (%d)"%len(x21)))
trace_x2 = go.Scatter(
       x = [x21.mean(),x21.mean()],
       y = [0,40],
       mode = "lines",
       name = "pool1_mean",
       )

fig.add_trace(trace_x2)

# Overlay both histograms
fig.update_layout(barmode='overlay')
# Reduce opacity to see both histograms
fig.update_traces(opacity=0.5)
fig.show()
