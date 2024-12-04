
import numpy as np
import networkx as nx
from NEMtropy import DirectedGraph
from NEMtropy import UndirectedGraph
from NEMtropy import matrix_generator as mg
from NEMtropy.network_functions import build_adjacency_from_edgelist
#from NEMtropy.network_functions import build_graph_from_edgelist
import random

# change working directory
import os
cwd = os.getcwd()
print(cwd)
os.chdir('/home/david/Work/Projects/NETMAP/data/unipartite_adjacency_matrices')
cwd = os.getcwd()
print(cwd)

list_unipartite_matrices_files = os.listdir()

list_unipartite_matrices_names = [x.replace('.csv', '') for x in list_unipartite_matrices_files]

# Number of random copies
number_random_samples = 100

random.seed(1235)

for mat_i in range(0, len(list_unipartite_matrices_files)):
    
    unip_mat_names = np.loadtxt('../unipartite_adjacency_matrices/' + list_unipartite_matrices_files[mat_i], delimiter=',',dtype=str)
    
    print(list_unipartite_matrices_files[mat_i])
    nodeto = unip_mat_names[1:,0]
    nodefrom = unip_mat_names[0, 1:]

    unip_mat = unip_mat_names[1:, 1:].astype(np.ubyte)
    nodeto_dict = dict(enumerate(np.unique(nodeto)))
    nodeto_inv_dict = {v:k for k,v in nodeto_dict.items()}
    nodefrom_dict = dict(enumerate(np.unique(nodefrom)))
    nodefrom_inv_dict = {v:k for k,v in nodefrom_dict.items()}

    edgelist = [(i,j) for i in range(unip_mat.shape[0]) for j in range(unip_mat.shape[0]) if unip_mat[i,j]!=0]
    edgelist_names = [(nodeto_dict[edge[0]], nodefrom_dict[edge[1]]) for edge in edgelist]

    dseq_in = unip_mat.sum(axis=0)
    dseq_out = unip_mat.sum(axis=1)

    graph = UndirectedGraph(unip_mat)
    # graph = DirectedGraph(unip_mat)

    # graph.solve_tool(model="dcm_exp",
    graph.solve_tool(model="cm_exp",
                    method="newton",
                    initial_guess="random")
    
    graph.ensemble_sampler(number_random_samples, cpu_n=5, output_dir="../null_unipartite_adjacency_matrices/")
    
    # Turning the edge lists into an adjacency matrices with node names 
    for sample_i in range(0, number_random_samples):
    
        # Load edge_list
        edgelist_ens = np.loadtxt("../null_unipartite_adjacency_matrices/"+repr(sample_i)+".txt")
        
        # Create an empty auxiliary matrix (zero matrix)
        new_unip_mat_names = unip_mat_names
        
        new_unip_mat_names[1:, 1:] = np.zeros(np.shape(new_unip_mat_names[1:, 1:]), dtype=int)
        
        # Fill in the empty auxiliary matrix. We consider that our matrices are directed
        for edge in edgelist_ens:
            nodo1, nodo2 = edge + 1
            new_unip_mat_names[int(nodo1)][int(nodo2)] = 1
            new_unip_mat_names[int(nodo2)][int(nodo1)] = 1 # For undirected matrices

        # Save the directed auxiliary matrix
        full_name = list_unipartite_matrices_files[mat_i] 
        network_name = full_name.replace(".csv", "")
        np.savetxt('../null_unipartite_adjacency_matrices/' + network_name + "_sample_" + repr(sample_i) + '.csv',new_unip_mat_names , fmt='%s', delimiter=',')
        
        # Remove the edge list
        os.remove("../null_unipartite_adjacency_matrices/"+repr(sample_i)+".txt")

