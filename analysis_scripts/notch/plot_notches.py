import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.figure()
df = df.from_csv("notch_elastic.csv")
plt.figure()
plt.plot(time,max_stress/1e6,label="stress 1")
plt.plot(time,damage,label="damage")
plt.legend()
plt.show()
