import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from paths import path

K_df = pd.read_csv(path('data', 'abundances', 'k_hypatia-28102024.csv'))
Eu_df = pd.read_csv(path('data', 'abundances', 'eu_hypatia-28102024.csv'))

df = K_df.merge(Eu_df, how = 'outer')

df['[K/Mg]'] = df['K'] - df['Mg']
df['[Eu/Mg]'] = df['Eu'] - df['Mg']
df['[K/Mg] err'] = df['K_err'] + df['Mg_err']
df['[Eu/Mg] err'] = df['Eu_err'] + df['Mg_err']
df['exohost'] = df['f_name'].notna()

# Filters out data points outside of the given percentile range.
def perc(data, min, max):
    return (data > np.percentile(data[data.notna()], min)) & (data < np.percentile(data[data.notna()], max))

fig = plt.figure(layout = 'constrained')
axis = fig.add_subplot()
fig.set_size_inches(5, 5)

exohosts = (df['exohost'] == True)
nonhosts = np.invert(exohosts)

K_perc = perc(df['[K/Mg]'], 1, 99) & df['[K/Mg]'].notna()
Eu_perc = perc(df['[Eu/Mg]'], 1, 99) & df['[Eu/Mg]'].notna()
mask = K_perc & Eu_perc

axis.errorbar(df['[K/Mg]'][nonhosts & mask], df['[Eu/Mg]'][nonhosts & mask], df['[K/Mg] err'][nonhosts & mask], df['[Eu/Mg] err'][nonhosts & mask], elinewidth = 0.5, alpha = 0.5, color = 'red', fmt = 'none', zorder = 1)
axis.scatter(df['[K/Mg]'][nonhosts & mask], df['[Eu/Mg]'][nonhosts & mask], s = 15, color = 'red', label = 'Non planet hosts')
axis.errorbar(df['[K/Mg]'][exohosts & mask], df['[Eu/Mg]'][exohosts & mask], df['[K/Mg] err'][exohosts & mask], df['[Eu/Mg] err'][exohosts & mask], elinewidth = 0.5, alpha = 0.5, color = 'blue', fmt = 'none', zorder = 1)
axis.scatter(df['[K/Mg]'][exohosts & mask], df['[Eu/Mg]'][exohosts & mask], s = 15, color = 'blue', label = 'Planet hosts')
axis.set_xlabel('[K/Mg] Relative Abundance [dex]', fontsize = 14)
axis.set_ylabel('[Eu/Mg] Relative Abundance [dex]', fontsize = 14)
axis.tick_params(labelsize = 12)
axis.legend(loc = 'upper left')

hist_ax_x = axis.inset_axes([0, 1.025, 1, 0.25], sharex = axis)
hist_ax_y = axis.inset_axes([1.025, 0, 0.25, 1], sharey = axis)
hist_ax_x.tick_params(axis = 'x', labelbottom = False)
hist_ax_y.tick_params(axis = 'y', labelleft = False)

hist_ax_x.hist(df['[K/Mg]'][K_perc][nonhosts], color = 'red', bins = 30)
hist_ax_x.set_ylabel('counts')
hist_ax_y.hist(df['[Eu/Mg]'][Eu_perc][nonhosts], orientation='horizontal', color = 'red', bins = 30)
hist_ax_y.set_xlabel('counts')
hist_ax_x.hist(df['[K/Mg]'][K_perc][exohosts], color = 'blue', bins = 15)
hist_ax_y.hist(df['[Eu/Mg]'][Eu_perc][exohosts], orientation='horizontal', color = 'blue', bins = 15)

fig.savefig(path('plots', 'eu_vs_k'), dpi = 300, bbox_inches = 'tight')