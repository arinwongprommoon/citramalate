import numpy as np
import pandas as pd
import gc

A = pd.DataFrame({'score':np.random.random(32346018), 'coordinate':np.arange(1, 32346019)})
peak_centers = []
peak_scores = []
exclusion = 147
count = 0
while A.shape[0] > 0:
    gc.collect()  # Force garbage collection.
    count += 1    # Increment the iteration count.
    print('iteration %d, shape %s' % (count, A.shape))
    raw_input()   # Wait for the user to press Enter.
    idx = A.score.values.argmax()
    one_center = A.coordinate.iloc[idx]
    # peak_centers and peak_scores are python lists
    peak_centers.append(one_center)
    peak_scores.append(A.score.iloc[idx])
    print(len(peak_centers), len(peak_scores))
    # exclude the coordinates around the selected peak
    A = A.loc[(A.coordinate <= one_center - exclusion) | (A.coordinate >= one_center + exclusion)]
