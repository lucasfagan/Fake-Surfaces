import numpy as np
import networkx as nx
from itertools import combinations, chain
# import pygraphviz as pgv
# from networkx.drawing.nx_agraph import graphviz_layout


SIZE_OF_MATRIX = 6

def find_lists(n, list_sum, partial_list=[]):
    if n == 1:
        # If we are at the last element, it needs to be equal to the remaining sum
        if list_sum >= 0:
            yield partial_list + [list_sum]
    elif n > 1:
        # Recursive step: try every possible value for the next element
        for i in range(list_sum + 1):
            yield from find_lists(n - 1, list_sum - i, partial_list + [i])

def create_adj_matrices(curr_matrix, all_matrices):
    # this function enumerates all symmetric matrices where the sum or each row (and thus column) is 4

    curr_row = len(curr_matrix)
    row_sum = 4
    for i in range(curr_row):
        row_sum -= curr_matrix[i][curr_row]
        if row_sum < 0:
            return
    # now we create all possible next rows
    rows = []

    # find all possible rows of 6-curr_row numbers that sum to row_sum
    # starting with [row_sum, 0, ...], [row_sum-1, 1, 0, ...], [row_sum-1, 0, 1, 0, ...] to [0, ..., row_sum]
    len_of_rows = SIZE_OF_MATRIX-curr_row
    for row in find_lists(len_of_rows, row_sum):
        rows.append(row)
    
    # prepend the items to make the row length 6 from the previous rows 
    beginning_of_row = []
    for i in range(curr_row):
        beginning_of_row.append(curr_matrix[i][curr_row])
    rows = [beginning_of_row + row for row in rows]


    # remove the rows that have a 4 or would cause a column sum of more than 4
    rows = [row for row in rows if ((max_sum_of_column(curr_matrix + [row]) <= 4) and (4 not in row) and (row[curr_row] % 2 == 0))]

    # recursive step
    for row in rows:
        if curr_row == SIZE_OF_MATRIX-1:
            all_matrices.append(curr_matrix + [row])
        else:
            create_adj_matrices(curr_matrix + [row], all_matrices)
    return all_matrices

    
def max_sum_of_column(matrix):
    # returns the maximum sum of any column of a matrix
    max_sum = 0
    for i in range(len(matrix)):
        curr_sum = 0
        for j in range(len(matrix)):
            curr_sum += matrix[j][i]
        if curr_sum > max_sum:
            max_sum = curr_sum
    
    return max_sum

def is_connected(matrix):
    n = len(matrix) # The number of nodes in the graph
    if n != SIZE_OF_MATRIX:
        raise ValueError("The matrix is not the right size")

    visited = [False] * n  # Array to keep track of visited nodes

    # The DFS function that marks the visited nodes
    def dfs(node):
        visited[node] = True
        for neighbor in range(n):
            if matrix[node][neighbor] != 0 and not visited[neighbor]:
                dfs(neighbor)

    # Perform DFS starting from the first node (index 0)
    dfs(0)

    # If DFS visited all nodes, the graph is connected
    return all(visited)

def unique_graphs_up_to_isomorphism(matrices):
    # Convert adjacency matrices to graph objects
    graphs = [nx.from_numpy_array(matrix, parallel_edges=True, create_using=nx.MultiGraph) for matrix in matrices]
    unique_graphs = []

    # Helper function to check if a graph is isomorphic to any in a list
    def is_new_isomorph(graph, graph_list):
        for g in graph_list:
            if nx.is_isomorphic(graph, g):
                return False
        return True

    # Loop over the graphs and check for isomorphisms
    i = 0
    for graph in graphs:

        i+=1
        if i % 10000 == 0:
            print(i)
        if is_new_isomorph(graph, unique_graphs):
            unique_graphs.append(graph)
    
    # Convert graphs back to adjacency matrices
    unique_matrices = [nx.to_numpy_array(g, dtype=int) for g in unique_graphs]
    # Return unique adjacency matrices
    return unique_matrices

def find_maximal_tree(adj_matrix):
    # Ensure the adjacency matrix is 6x6
    if len(adj_matrix) != 6 or any(len(row) != 6 for row in adj_matrix):
        raise ValueError("Input is not a 6x6 adjacency matrix")
    
    # Initialize variables
    num_vertices = 6
    visited = [False] * num_vertices
    tree_edges = []
    
    # Recursive DFS function to find the maximal tree
    def dfs(current_vertex, parent_vertex):
        visited[current_vertex] = True
        for neighbor in range(num_vertices):
            if adj_matrix[current_vertex][neighbor] > 0:  # Check for an edge
                # Decrease the edge count in the adjacency matrix
                adj_matrix[current_vertex][neighbor] -= 1
                adj_matrix[neighbor][current_vertex] -= 1
                # If the neighbor has not been visited, it's a tree edge
                if not visited[neighbor]:
                    tree_edges.append((current_vertex + 1, neighbor + 1))  # +1 for 1-indexed vertices
                    dfs(neighbor, current_vertex)
                # If the neighbor has been visited and is not the parent, we're encountering a back edge, ignore to prevent cycles
    
    # Starting DFS from vertex 0 (1-indexed as vertex 1)
    dfs(0, -1)
    
    # Check if all vertices were included in the tree, making it a maximal tree
    if len(tree_edges) != num_vertices - 1:
        raise ValueError("The input graph does not have a spanning tree.")

    return tree_edges

def generate_graph_by_edges(adj_matrix):
    graph_by_edges = {'x':[], 'y':[], 'z':[],'u':[],'v':[],'w':[]}
    # loop through the upper right triangle of the adjaency matrix counting diagonals (i.e., the edges). 
    # each edge needs to be given a number from 1 to 12, and added to its in and out vertex
    # for example, if the first nonzero entry is a 2 in the first row, third column, we add 1,2 to the x list, and -1,-2 to the z list
    curr_edge_to_add = 1
    verts_order = ['x','y','z','u','v','w']
    for i in range(len(adj_matrix)):
        for j in range(i,len(adj_matrix)):
            if i==j and adj_matrix[i][j] == 2:
                # add the edge to the appropriate lists
                graph_by_edges[verts_order[i]].append(curr_edge_to_add)
                graph_by_edges[verts_order[j]].append(-1*curr_edge_to_add)
                curr_edge_to_add += 1
            else:
                for k in range(adj_matrix[i][j]):
                    # add the edge to the appropriate lists
                    graph_by_edges[verts_order[i]].append(curr_edge_to_add)
                    graph_by_edges[verts_order[j]].append(-1*curr_edge_to_add)
                    curr_edge_to_add += 1

    for key in graph_by_edges.keys():
        graph_by_edges[key] = tuple(graph_by_edges[key])

    return graph_by_edges

def change_format_maxl_tree(tree_edges, graph_by_edges):
    # right now, tree edges contains things like (1,4). We want to change it to a single number, representing the edge number
    # so in the specific example of (1,4), we want to find an edge "i" such that i is in the x list and -i is in the u list

    # first, change the tree_edges to be 0-indexed
    tree_edges = [(edge[0]-1, edge[1]-1) for edge in tree_edges]

    verts_order = ['x','y','z','u','v','w']
    # loop through the edges, and find the edge number
    tree_edges_numbers = []
    for edge in tree_edges:
        # find the edge number
        vertex1 = verts_order[edge[0]]
        vertex2 = verts_order[edge[1]]
        # find intersection of the absolute values of the two lists
        edges1 = [abs(edge) for edge in graph_by_edges[vertex1]]
        edges2 = [abs(edge) for edge in graph_by_edges[vertex2]]
        edge_num = list(set(edges1).intersection(edges2))[0]
        tree_edges_numbers.append(edge_num)
    return tree_edges_numbers


def generate_by_verts(graph_by_edges):
    graph_by_verts = {}
    # combine all the values for all the keys in graph_by_edges into a list, and get the max
    num_edges = max([abs(edge) for edge in list(chain(*graph_by_edges.values()))])
    for i in range(1,num_edges+1):
        # have i map to whichever key maps to i in graph_by_edges
        i_keys = [key for key in graph_by_edges.keys() if i in graph_by_edges[key]]
        if len(i_keys) != 1:
            print("graph_by_edges:",graph_by_edges)
            print("i:",i)
            exit()
        graph_by_verts[i] = i_keys[0]
        minus_i_keys = [key for key in graph_by_edges.keys() if -1*i in graph_by_edges[key]]
        if len(minus_i_keys) != 1:
            print("graph_by_edges:",graph_by_edges)
            print("-i:",-1*i)
            exit()
        graph_by_verts[-1*i] = minus_i_keys[0]
    return graph_by_verts


all_matrices = create_adj_matrices([], [])
# only keep connected matrices, and turn them into numpy arrays
all_connected_matrices_np = [np.array(matrix) for matrix in all_matrices if is_connected(matrix)]
print(len(all_connected_matrices_np))
all_unique_matrices_np = unique_graphs_up_to_isomorphism(all_connected_matrices_np)
print(all_unique_matrices_np[:5])
print(len(all_unique_matrices_np))

for matrix in all_unique_matrices_np:
    graph_by_edges = generate_graph_by_edges(matrix)
    graph_by_verts = generate_by_verts(graph_by_edges)
    maxl_tree = change_format_maxl_tree(find_maximal_tree(matrix), graph_by_edges)
    print("graph_six_verts_"+str(matrix_num)+" = (",end="")
    print(graph_by_verts, end="")
    print(",",end="")
    print(graph_by_edges, end="")
    print(",",end="")
    print(maxl_tree, end="")
    print(")")

   


    
