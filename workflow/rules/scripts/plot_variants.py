#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt


with open(sys.argv[1], "r") as file_in:
    # Load the data into a pandas dataframe
    df = pd.read_table('file_in')


# Plot each feature in a separate subplot
fig, ax = plt.subplots(3, 2, figsize=(15, 10))
df['QD'].plot.hist(ax=ax[0][0], bins=100)
ax[0][0].set_title('QD')

df['FS'].plot.hist(ax=ax[0][1], bins=100)
ax[0][1].set_title('FS')

df['SOR'].plot.hist(ax=ax[1][0], bins=100)
ax[1][0].set_title('SOR')

df['MQ'].plot.hist(ax=ax[1][1], bins=100)
ax[1][1].set_title('MQ')

df['MQRankSum'].plot.hist(ax=ax[2][0], bins=100)
ax[2][0].set_title('MQRankSum')

df['ReadPosRankSum'].plot.hist(ax=ax[2][1], bins=100)
ax[2][1].set_title('ReadPosRankSum')

plt.tight_layout()
plt.show()
