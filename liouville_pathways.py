import numpy as np
import itertools

n_states = 4 #including initial (e.g. ground) state
n_lasers = n_states-1
t_orderings = np.array(list(itertools.permutations(np.arange(n_lasers), n_lasers)))
pm = np.array([1, -1, 1])
kb = np.array([1, -1])
print(kb.reshape(1,len(kb)).T*pm.reshape(1, len(pm))*t_orderings)

state_labels = ['0', '1', '2', '3']


        
        






