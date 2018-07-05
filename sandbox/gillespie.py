N = 350 # total population
T = 100.0 # maximum elapsed time
t = 0.0 # start time
V = 100.0 # spatial parameter
alpha = 10.0 # rate of infection after contact
beta = 0.5 # rate of cure
n_I = 1 # initial infected population

import numpy as np
# set numercial types to numpy words
N = np.int64(N)
T = np.float(T)
t = np.float(t)
V = np.float64(V)
alpha = np.float64(alpha)
beta = np.float64(beta)
n_I = np.int64(n_I)

# compute susceptible population, set recovered to zero
n_S = N - n_I
n_R = np.int64(0)

# initialise results list
SIR_data = []
SIR_data.append((t, n_S, n_I, n_R))

#main loop
while t < T:

    if n_I == 0:
        break

    w1 = alpha * n_S * n_I / V
    w2 = beta * n_I
    W = w1 + w2

    dt = np.log(np.random.random_sample()) / W
    t = t + dt

    if np.random.random_sample() < w1 / W:
        n_S = n_S - 1
        n_I = n_I + 1

    else:
        n_I = n_I - 1
        n_R = n_R + 1

    SIR_data.append((t, n_S, n_I, n_R))

with open('SIR_data.txt', 'w+') as fp:
    fp.write('\n'.join('%f %i %i %i' % x for x in SIR_data))
