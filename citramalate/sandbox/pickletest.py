import cPickle as pickle
import numpy as np

data = [np.arange(8).reshape(2, 4), np.arange(10).reshape(2, 5)]

with open('mat.pkl', 'wb') as outfile:
    pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)

with open('mat.pkl', 'rb') as infile:
    result = pickle.load(infile)
