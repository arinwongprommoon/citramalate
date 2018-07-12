import numpy as np
import pandas as pd

A = pd.DataFrame({'score':np.random.random(132346018), 'coordinate':np.arange(1, 132346019)})
peak_centers = []
peak_scores = []
exclusion = 147
while A.shape[0] > 0:
    idx = A.score.values.argmax()
    one_center = A.coordinate.iloc[idx]
    # peak_centers and peak_scores are python lists
    peak_centers.append(one_center)
    peak_scores.append(A.score.iloc[idx])
    # exclude the coordinates around the selected peak
    A = A.loc[(A.coordinate <= one_center - exclusion) | (A.coordinate >= one_center + exclusion)]

# terminated the loop after memory consumption gets to 90% of machine memory
# but peak_centers and peak_scores are still short lists
print len(peak_centers)
# output is 16
