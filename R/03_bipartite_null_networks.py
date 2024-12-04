
import sys
sys.path.append('..')
from bicm import *
import numpy as np
import random
import contextlib

# change working directory
import os
cwd = os.getcwd()
print(cwd)
os.chdir('/home/david/Work/Projects/NETMAP/data/bipartite_adjacency_matrices')
cwd = os.getcwd()
print(cwd)


list_bipartite_matrices_files = os.listdir()


list_bipartite_matrices_names = [x.replace('.csv', '') for x in list_bipartite_matrices_files]

# Number of random copies
number_random_samples = 100


random.seed(123)

for mat_i in range(0, len(list_bipartite_matrices_files)):
    
    biad_mat_names = np.loadtxt('../bipartite_adjacency_matrices/' + list_bipartite_matrices_files[mat_i], delimiter=',',dtype=str)
    #print(list_bipartite_matrices_files[mat_i])
    plants = biad_mat_names[1:,0]
    pollinators = biad_mat_names[0, 1:]

    biad_mat = biad_mat_names[1:, 1:].astype(np.ubyte)
    plants_dict = dict(enumerate(np.unique(plants)))
    plants_inv_dict = {v:k for k,v in plants_dict.items()}
    pollinators_dict = dict(enumerate(np.unique(pollinators)))
    pollinators_inv_dict = {v:k for k,v in pollinators_dict.items()}


    edgelist = edgelist_from_biadjacency(biad_mat)[0]
    edgelist_names = [(plants_dict[edge[0]], pollinators_dict[edge[1]]) for edge in edgelist]
    
    try:
        myGraph = BipartiteGraph(edgelist=edgelist_names)
        myGraph = BipartiteGraph()
        myGraph.set_edgelist(edgelist_names)
        with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
          myGraph.solve_tool(verbose=False)
        
        dict_x, dict_y = myGraph.get_bicm_fitnesses()
        #print('Yielded data type is:', type(dict_x))

        #myGraph.loglikelihood

        x = myGraph.x
        y = myGraph.y
        rows_dict = myGraph.rows_dict
        cols_dict = myGraph.cols_dict
        avg_mat = myGraph.get_bicm_matrix()
        #print(avg_mat[0, 0] == x[0] * y[0] / (1 + x[0] * y[0]))

        avg_mat = bicm_from_fitnesses(x, y)
        #print(avg_mat[0, 0] == x[0] * y[0] / (1 + x[0] * y[0]))

        #myGraph.check_sol(biad_mat, avg_mat)

        sample_with_names = np.zeros_like(biad_mat_names)
        sample_with_names[:, 0] = biad_mat_names[:,0]
        sample_with_names[0, :] = biad_mat_names[0,:]
        for sample_i in range(number_random_samples):

            sampled_network = sample_bicm(avg_mat)

            while True:
                try:
                    sample_with_names[1:,1:] = sampled_network
                except ValueError:
                    continue
                break
            full_name = list_bipartite_matrices_files[mat_i] 
            network_name = full_name.replace(".csv", "")
            np.savetxt('../null_bipartite_adjacency_matrices/' + network_name + "_sample_" + repr(sample_i) + '.csv', sample_with_names, fmt='%s', delimiter=',')
    
    except IndexError as e:
        print("************ error : "+list_bipartite_matrices_files[mat_i])

