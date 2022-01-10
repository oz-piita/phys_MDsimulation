import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("./result.csv",header=None)
df.columns = ["time","kinetic energy","potenstial energy","total energy"]

df = df.drop("time",axis=1)

plt.figure()
df.plot()
plt.show()