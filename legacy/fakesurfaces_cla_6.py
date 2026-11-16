from itertools import combinations, product, permutations, chain
from numpy.linalg import det
from numpy import matrix
from collections import Counter
import time
import sys

def sign(x):
    return (x > 0) - (x < 0)

def dfs(sequence, valid_sequences=None, graph=None):
    graph_by_edges, graph_by_vertices, maxl_tree = graph
    if valid_sequences is None:    #Set default value inside function
        valid_sequences = set()
    if len(sequence)>1 and sequence[-1]==sequence[0]:
        if len(sequence)>3 and len(sequence)<20:
        # valid_sequences.append(' -> '.join(str(s) for s in sequence))
        # valid_sequences.append(list(sequence)[:-1]))
            valid_sequences.add(tuple(disk_to_std_form(list(sequence)[:-1])))
            return valid_sequences
    incoming_edge, vertex, outgoing_edge = sequence[-1] # this is the last item in the sequence
    incoming_edge_new = -1*outgoing_edge #since we are on the other side of the edge
    destination_vertex = graph_by_edges[incoming_edge_new] 
    possible_next_edges = [edge for edge in list(graph_by_vertices[destination_vertex]) if edge != incoming_edge_new]
    for edge in possible_next_edges:
        if (incoming_edge_new, destination_vertex, edge) in sequence[1:] or (edge, destination_vertex, incoming_edge_new) in sequence[1:]:
            # edge pair already been used for that vertex, so
            continue
        sequence.append((incoming_edge_new, destination_vertex, edge))
        result = dfs(sequence, valid_sequences, graph) # Save returned value in a variable
        if result:   # Check if it is not None
            valid_sequences = result
        # backtrack
        sequence.pop()
    return valid_sequences  # Return at the end of function

def disk_to_std_form(disk, fwd_or_backwards=1):
    min_val_idx = min(enumerate(disk), key=lambda x: x[1])[0]
    if fwd_or_backwards:
        return disk[min_val_idx:] + disk[:min_val_idx]
    else:
        return disk[min_val_idx::-1] + disk[:min_val_idx:-1]

def find_partitions(sequence_lengths, num_disks, k):
    result = []
    tried_combos = []
    for combination in product(*([sequence_lengths]*num_disks)):
        if sum(combination)==k:
            if sorted(combination) in tried_combos:
                continue
            else:
                tried_combos.append(sorted(combination))
                result.append(combination)

    return result

def create_disks(graph):
    graph_by_edges, graph_by_vertices, maxl_tree = graph
    disks = set()
    for vertex in graph_by_vertices.keys():
        for edge1,edge2 in combinations(graph_by_vertices[vertex],2):
            print("Trying edges",edge1,edge2)
            # add on the new one
            disks = dfs([(edge1,vertex,edge2)], disks, graph)
    return disks

def remove_duplicate_disks(disks):
    unique_disks = set()
    for disk in disks:
        disk = tuple(disk_to_std_form(disk))
        disk_backwards = tuple(disk_to_std_form([(x[2],x[1],x[0]) for x in disk][::-1]))
        if disk not in unique_disks and disk_backwards not in unique_disks:
            unique_disks.add(disk)
    return unique_disks

def split_up_sequences_by_length_and_sort(sequences):
    sequence_lengths = {}
    sequences_sorted = {}
    for sequence in sequences:
        sequences_sorted[sequence] = [tuple(sorted((elt[0],elt[2]))) for elt in sequence]
        if len(sequence) not in sequence_lengths.keys():
            sequence_lengths[len(sequence)] = [sequence]
        else:
            sequence_lengths[len(sequence)].append(sequence)
    return (sequence_lengths, sequences_sorted)

def find_det_of_B(fake_surface, graph):
    graph_by_edges, graph_by_vertices, maxl_tree = graph
    B = []
    for disk in fake_surface:
        coefficent_dict = {}
        for edge in [x[0] for x in disk]: 
            if abs(edge) not in coefficent_dict.keys():
                coefficent_dict[abs(edge)] = sign(edge)
            else:
                coefficent_dict[abs(edge)] += sign(edge)
        for edge in sorted(list(set([abs(x) for x in graph_by_edges.keys()]))):
            if edge not in coefficent_dict.keys():
                coefficent_dict[edge] = 0
        # collapse maximal tree
        for edge in maxl_tree:
            del coefficent_dict[edge]
        disks_column_in_B_matrix = [coefficent_dict[edge] for edge in sorted(coefficent_dict.keys())]
        B.append(disks_column_in_B_matrix)
    return round(det(matrix(B)))

def find_acyclic_fake_surfaces(sequences, graph):
    graph_by_edges, graph_by_vertices, maxl_tree = graph
    sequence_lengths, sequences_sorted = split_up_sequences_by_length_and_sort(sequences)
    # print("Length of unique sequences",{key:len(sequence_lengths[key]) for key in sequence_lengths.keys()})
    acyclic_fake_surfaces = set()
    all_partitions = set()
    acyclic_partitions = set()
    num_verts = len(graph_by_vertices.keys())
    num_disks, total_edge_pairs = num_verts+1, 6*num_verts
    for partition in find_partitions(sequence_lengths.keys(), num_disks, total_edge_pairs):
        print("Trying partition",partition)
        sequences_of_right_lengths = [sequence_lengths[length] for length in partition]
        for choice_of_seqs in product(*sequences_of_right_lengths):
            union = set()
            for sequence in choice_of_seqs:
                union.update(sequences_sorted[sequence])
            # need to sort so (-4,1) is same as (4,-1)
            if len((union)) == total_edge_pairs: # we already know it's the right length
                all_partitions.add(tuple(partition))
                if abs(find_det_of_B(choice_of_seqs, graph)) == 1:
                    acyclic_partitions.add(tuple(partition))
                    acyclic_fake_surfaces.add(tuple(sort_surface_longest_disk_to_shortest(choice_of_seqs)))
    return (acyclic_fake_surfaces, all_partitions, acyclic_partitions)

def write_to_file_for_sage(fake_surfaces,graph,graph_num):
    graph_by_edges, graph_by_vertices, maxl_tree = graph
    # fake_surfaces = [[[edge[0] for edge in disk] for disk in surface] for surface in fake_surfaces]
    #TODO: MAKE SURE THIS FUNCTION WORKS, IN PARTICULAR, THE LINE ABOVE SHOULD ACTUALLY BE COMMENTED OUT 
    fake_surfaces = [[[elt for elt in disk if (elt not in maxl_tree and -1*elt not in maxl_tree)] for disk in surface] for surface in fake_surfaces]
    # sort maxl_tree largest to smallest
    for removed_edge in sorted(maxl_tree, reverse=True):
        fake_surfaces = [[[edge-1 if edge>removed_edge else edge for edge in disk] for disk in surface] for surface in fake_surfaces]
        fake_surfaces = [[[edge+1 if edge<-1*removed_edge else edge for edge in disk] for disk in surface] for surface in fake_surfaces]
    with open('/home/lucasfagan/fakesurfaces_complexity6/groups_graph'+str(graph_num)+'.txt', 'w') as f:
        for surface in fake_surfaces:
            f.write(str(surface)+'\n')

def sort_surface_longest_disk_to_shortest(t):
    return sorted(sorted(t)[::-1], key=lambda x: len(x), reverse=True)

def remove_duplicate_fake_surfaces(fake_surfaces):
    unique_fake_surfaces = {}
    for fake_surface in fake_surfaces:
        # fake_surface = tuple(sort_surface_longest_disk_to_shortest([disk_to_std_form(tuple([edge[0] for edge in disk])) for disk in fake_surface]))
        fake_surface = tuple(sort_surface_longest_disk_to_shortest([(tuple([edge[0] for edge in disk])) for disk in fake_surface]))
        partition_used = tuple(sorted([len(sequence) for sequence in fake_surface]))
        if partition_used not in unique_fake_surfaces.keys():
            unique_fake_surfaces[partition_used] = set([fake_surface])
        else:
            already_in = False
            for already_seen in unique_fake_surfaces[partition_used]:
                # if are_same(fake_surface, already_seen) and not are_same(already_seen, fake_surface):
                #     print("Found a duplicate surface that is not the same!")
                #     print("Surface 1:",fake_surface)
                #     print(are_same(fake_surface, already_seen))
                #     print("Surface 2:",already_seen)
                if are_same(fake_surface, already_seen):
                    already_in = True
                    break
            if not already_in:
                unique_fake_surfaces[partition_used].add(fake_surface)
    return unique_fake_surfaces

def are_compatible(disk1, disk2, d):
    starting_d = d.copy()
    valid_dicts = []
    l = len(disk1)
    for shift_amount in range(l):
        for sgn in [-1,1]:
            d = starting_d.copy()
            for i in range(len(disk1)):
                # check: does a key map to the wrong value or is a value already mapped to by a different key?
                if disk1[i] in d.keys() and d[disk1[i]] != disk2[(shift_amount+sgn*i) % len(disk1)]:
                    break
                if disk2[(shift_amount+sgn*i) % len(disk1)] in d.values() and disk1[i] != list(d.keys())[list(d.values()).index(disk2[(shift_amount+sgn*i) % len(disk1)])]:
                    break
                d[disk1[i]] = disk2[(shift_amount+i*sgn) % len(disk1)]
            else:
                valid_dicts.append(d)
    return valid_dicts

def reorder_disks_of_same_size(surface):
    len_dict = {}
    for sublist in surface:
        len_dict.setdefault(len(sublist), []).append(sublist)
    permutations_list = [list(permutations(len_dict[key])) for key in len_dict.keys()]
    final_list = [tuple(chain(*item)) for item in product(*permutations_list)]
    return final_list

def are_same(surface1, surface2):
    for surface2_version in reorder_disks_of_same_size(surface2):
        if are_same_with_ordering(surface1,surface2_version):
            return True
    return False

def are_same_with_ordering(surface1,surface2):
    for disk_num, disk in enumerate(surface1):
        if len(Counter(disk)) != len(Counter(surface2[disk_num])):
            return False
    # now try to transform disk1 into disk2
    transform_dicts = [{}]
    for disk_num in range(len(surface1)): #iterate through disks
        new_transform_dicts = []
        for transform_dict in transform_dicts:
            disk1 = surface1[disk_num]
            disk1_neg = tuple([-1*x for x in disk1])
            disk2 = surface2[disk_num]
            disk2_neg = tuple([-1*x for x in disk2])
            temp_dicts = are_compatible(disk1, disk2, transform_dict)
            if temp_dicts:
                new_transform_dicts.extend(temp_dicts)
            temp_dicts = are_compatible(disk1, disk2_neg, transform_dict)
            if temp_dicts:
                new_transform_dicts.extend(temp_dicts)
            temp_dicts = are_compatible(disk1_neg, disk2, transform_dict)
            if temp_dicts:
                new_transform_dicts.extend(temp_dicts)
            temp_dicts = are_compatible(disk1_neg, disk2_neg, transform_dict)
            if temp_dicts:
                new_transform_dicts.extend(temp_dicts)
        transform_dicts = new_transform_dicts
    if len(transform_dicts) == 0:
        return False
    return True

graph_four_verts_1 = ({1:'x',-1:'z',2:'x',-2:'z',3:'x',-3:'y',4:'x',-4:'y',5:'y',-5:'w',6:'y',-6:'w',7:'w',-7:'z',8:'w',-8:'z'},{'x':(1,2,3,4),'y':(-3,-4,5,6),'z':(-1,-2,-7,-8),'w':(-5,-6,7,8)},[1,3,5]) # 8 hours
graph_four_verts_2 = ({1:'x',-1:'x',2:'y',-2:'y',3:'z',-3:'z',4:'w',-4:'w',5:'x',-5:'y',6:'y',-6:'z',7:'z',-7:'w',8:'x',-8:'w'},{'x':(1,-1,5,8),'y':(2,-2,-5,6),'z':(3,-3,-6,7),'w':(4,-4,-7,-8)},[5,6,7]) # fastish
graph_four_verts_3 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'y',4:'y',-4:'z',5:'y',-5:'z',6:'z',-6:'w',7:'z',-7:'w',8:'w',-8:'w'},{'x':(1,-1,2,3),'y':(-2,-3,4,5),'z':(-4,-5,6,7),'w':(-6,-7,8,-8)},[2,4,6]) # fastish
graph_four_verts_4 = ({1:'x',-1:'x',2:'y',-2:'y',3:'x',-3:'y',4:'x',-4:'z',5:'y',-5:'w',6:'z',-6:'w',7:'z',-7:'w',8:'z',-8:'w'},{'x':(1,-1,3,4),'y':(2,-2,-3,5),'z':(-4,6,7,8),'w':(-5,-6,-7,-8)},[3,5,6]) # ~2 hours
graph_four_verts_5 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'z',4:'y',-4:'y',5:'y',-5:'z',6:'z',-6:'w',7:'z',-7:'w',8:'w',-8:'w'},{'x':(1,-1,2,3),'y':(-2,4,-4,5),'z':(-3,-5,6,7),'w':(-6,-7,8,-8)},[2,5,6]) # ~2 hours
graph_four_verts_6 = ({1:'x',-1:'x',2:'x',-2:'y',3:'y',-3:'z',4:'y',-4:'w',5:'y',-5:'w',6:'x',-6:'w',7:'z',-7:'z',8:'z',-8:'w'},{'x':(1,-1,2,6),'y':(-2,3,4,5),'z':(-3,7,-7,8),'w':(-5,-4,-6,-8)},[2,3,8]) # ~1.5 hr
graph_four_verts_7 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'w',4:'y',-4:'z',5:'y',-5:'z',6:'y',-6:'w',7:'z',-7:'w',8:'z',-8:'w'},{'x':(1,-1,2,3),'y':(-2,4,5,6),'z':(-4,-5,7,8),'w':(-3,-6,-7,-8)},[2,4,7]) # ~1.5 hr
graph_four_verts_8 = ({1:'x',-1:'y',2:'x',-2:'z',3:'x',-3:'z',4:'x',-4:'w',5:'y',-5:'z',6:'y',-6:'w',7:'y',-7:'w',8:'z',-8:'w'},{'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-2,-3,-5,8),'w':(-4,-6,-7,-8)},[1,5,8]) # 1.5-2 hr
graph_four_verts_9 = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'y',4:'x',-4:'w',5:'y',-5:'z',6:'z',-6:'w',7:'z',-7:'w',8:'z',-8:'w'},{'x':(1,2,3,4),'y':(-1,-2,-3,5),'z':(-5,6,7,8),'w':(-4,-6,-7,-8)},[1,5,6]) # super long 
graph_four_verts_10 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'y',4:'y',-4:'z',5:'y',-5:'w',6:'z',-6:'w',7:'z',-7:'w',8:'z',-8:'w'},{'x':(1,-1,2,3),'y':(-2,-3,4,5),'z':(-4,6,7,8),'w':(-5,-6,-7,-8)},[2,5,6]) #1.5 hr
# graphs with 3 vertices
graph_three_verts_1 = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'z',4:'x',-4:'z',5:'y',-5:'z',6:'y',-6:'z'},{'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-3,-4,-5,-6)},[1,3])
graph_three_verts_2 = ({1:'x',-1:'x',2:'y',-2:'y',3:'z',-3:'z',4:'x',-4:'y',5:'y',-5:'z',6:'x',-6:'z'},{'x':(1,-1,4,6),'y':(2,-2,-4,5),'z':(3,-3,-5,-6)},[4,5])
graph_three_verts_3 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'y',4:'y',-4:'z',5:'y',-5:'z',6:'z',-6:'z'},{'x':(1,-1,2,3),'y':(-2,-3,4,5),'z':(-4,-5,6,-6)},[2,4])
graph_three_verts_4 = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'y',4:'x',-4:'z',5:'y',-5:'z',6:'z',-6:'z'},{'x':(1,2,3,4),'y':(-1,-2,-3,5),'z':(-4,-5,6,-6)},[1,4])
# graphs with two vertices:
graph_1_paper = ({1:'x',-1:'x',2:'y',-2:'y',3:'x',-3:'y',4:'x',-4:'y'},{'x':(-1,1,3,4),'y':(2,-2,-3,-4)},[3])
graph_2_paper = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'y',4:'x',-4:'y'},{'x':(1,2,3,4),'y':(-1,-2,-3,-4)},[1])
# graphs with one vertex:
graph_one_vert = ({1:'x',-1:'x',2:'x',-2:'x'},{'x':(1,-1,2,-2)},[])
# graphs with five vertices:
# graph_five_verts_1 = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'w',4:'x',-4:'w',5:'y',-5:'z',6:'y',-6:'z',7:'z',-7:'v',8:'z',-8:'v',9:'v',-9:'w',10:'v',-10:'w'},{'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-5,-6,7,8),'v':(-7,-8,9,10),'w':(-3,-4,-9,-10)},[1,5,8,9]) #double pentagon 
# graph_five_verts_2 =  ({1:'x',-1:'y',2:'x',-2:'z',3:'x',-3:'w',4:'x',-4:'v',5:'y',-5:'v',6:'y',-6:'w',7:'y',-7:'z',8:'z',-8:'v',9:'z',-9:'w',10:'w',-10:'v'},{'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-2,-7,8,9),'v':(-4,-5,-8,-10),'w':(-3,-6,-9,10)},[1,7,9,10]) #complete graph on 5 vertices
# graph_five_verts_3 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'y',4:'y',-4:'z',5:'y',-5:'z',6:'z',-6:'w',7:'z',-7:'w',8:'w',-8:'v',9:'w',-9:'v',10:'v',-10:'v'},{'x':(1,-1,2,3),'y':(-2,-3,4,5),'z':(-4,-5,6,7),'w':(-6,-7,8,9),'v':(-8,-9,10,-10)},[2,4,6,8]) #long chain graph
# graph_five_verts_4 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'v',4:'y',-4:'z',5:'y',-5:'z',6:'z',-6:'z',7:'v',-7:'w',8:'v',-8:'w',9:'w',-9:'w',10:'y',-10:'v'},{'x':(1,-1,2,3),'y':(-2,4,5,10),'z':(-4,-5,6,-6),'v':(-3,7,8,-10),'w':(-7,-8,9,-9)},[4,2,3,7]) # triangle with tails 
# graph_five_verts_5 = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'y',4:'x',-4:'z',5:'y',-5:'z',6:'z',-6:'v',7:'z',-7:'v',8:'v',-8:'w',9:'v',-9:'w',10:'w',-10:'w'},{'x':(1,2,3,4),'y':(-1,-2,-3,5),'z':(-4,-5,6,7),'v':(-6,-7,8,9),'w':(-8,-9,10,-10)},[1,5,6,8]) # pentagon with tails
# graph_five_verts_6 = ({1:'x',-1:'y',2:'x',-2:'y',3:'x',-3:'v',4:'x',-4:'w',5:'y',-5:'z',6:'y',-6:'z',7:'z',-7:'w',8:'z',-8:'v',9:'v',-9:'w',10:'v',-10:'w'},{'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-5,-6,7,8),'v':(-3,-8,9,10),'w':(-4,-7,-9,-10)},[1,5,8,9])  
# graph_five_verts_7 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'z',4:'y',-4:'y',5:'y',-5:'z',6:'z',-6:'w',7:'z',-7:'w',8:'w',-8:'v',9:'w',-9:'v',10:'v',-10:'v'},{'x':(1,-1,2,3),'y':(-2,4,-4,5),'z':(-3,-5,6,7),'v':(-8,-9,10,-10),'w':(-6,-7,8,9)},[2,5,6,8])  
# graph_five_verts_8 = ({1:'x',-1:'x',2:'x',-2:'y',3:'x',-3:'z',4:'y',-4:'w',5:'y',-5:'w',6:'z',-6:'z',7:'z',-7:'w',8:'y',-8:'v',9:'w',-9:'v',10:'v',-10:'v'},{'x':(1,-1,2,3),'y':(-2,4,5,8),'z':(-3,6,-6,7),'v':(-8,-9,10,-10),'w':(-4,-5,-7,9)},[2,3,7,9])  
# graph_five_verts_9 = ({1:'x',-1:'y',2:'x',-2:'z',3:'x',-3:'z',4:'x',-4:'w',5:'y',-5:'z',6:'y',-6:'w',7:'y',-7:'v',8:'z',-8:'w',9:'w',-9:'v',10:'v',-10:'v'},{'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-2,-3,-5,8),'v':(-7,-9,10,-10),'w':(-4,-6,-8,9)},[1,2,8,9])  

# graph = graph_one_vert

# all_disks = create_disks(graph) #dfs part 
# unique_disks = remove_duplicate_disks(all_disks)
# acyclic_fake_surfaces, all_partitions, acyclic_partitions = find_acyclic_fake_surfaces(unique_disks, graph)
# unique_acyclic_fake_surfaces = remove_duplicate_fake_surfaces(acyclic_fake_surfaces)

# # print(unique_acyclic_fake_surfaces)
# # print("Number of unique acyclic fake surfaces:",sum([len(unique_acyclic_fake_surfaces[key]) for key in unique_acyclic_fake_surfaces.keys()]))
# # print("Unique acyclic fake surface partitions:",{key:len(unique_acyclic_fake_surfaces[key]) for key in unique_acyclic_fake_surfaces.keys()})
# with open('/home/lucasfagan/unique_surfaces_graph_five_verts_9.txt', 'w') as f:
#     f.write("Number of unique acyclic fake surfaces: "+str(sum([len(unique_acyclic_fake_surfaces[key]) for key in unique_acyclic_fake_surfaces.keys()]))+'\n')
#     f.write(str(unique_acyclic_fake_surfaces)+'\n')

# fake_surfaces_no_dict = []
# for key in unique_acyclic_fake_surfaces.keys():
#     fake_surfaces_no_dict.extend(unique_acyclic_fake_surfaces[key])
# write_to_file_for_sage(fake_surfaces_no_dict, graph)

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

graph_five_verts_1_by_edges = {'x':(1,-1,2,3),'y':(-2,-3,4,5),'z':(-4,-5,6,7),'v':(-6,-7,8,9),'w':(-8,-9,10,-10)}
graph_five_verts_1 = (generate_by_verts(graph_five_verts_1_by_edges),graph_five_verts_1_by_edges,[2,4,6,8])
graph_five_verts_2_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-3,-5,-6,7),'v':(-4,-7,8,9),'w':(-8,-9,10,-10)}
graph_five_verts_2 = (generate_by_verts(graph_five_verts_2_by_edges),graph_five_verts_2_by_edges,[1,5,7,8])
graph_five_verts_3_by_edges = {'x':(1,-1,2,3),'y':(-2,4,5,6),'z':(-4,7,-7,8),'v':(-5,-8,9,10),'w':(-3,-6,-9,-10)}
graph_five_verts_3 = (generate_by_verts(graph_five_verts_3_by_edges),graph_five_verts_3_by_edges,[2,4,8,9])
graph_five_verts_4_by_edges = {'x':(1,2,3,4),'y':(-1,5,6,-6),'z':(-2,-5,7,8),'v':(-3,-7,9,10),'w':(-4,-8,-9,-10)}
graph_five_verts_4 = (generate_by_verts(graph_five_verts_4_by_edges),graph_five_verts_4_by_edges,[1,5,7,9])
graph_five_verts_5_by_edges = {'x':(1,2,3,4),'y':(-1,5,6,-6),'z':(-2,-5,7,-7),'v':(-8,-9,10,-10),'w':(-3,-4,8,9)}
graph_five_verts_5 = (generate_by_verts(graph_five_verts_5_by_edges),graph_five_verts_5_by_edges,[1,5,3,8])
graph_five_verts_6_by_edges = {'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-2,-5,-6,-7),'v':(-8,-9,10,-10),'w':(-3,-4,8,9)}
graph_five_verts_6 = (generate_by_verts(graph_five_verts_6_by_edges),graph_five_verts_6_by_edges,[1,5,3,8])
graph_five_verts_7_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-2,-6,7,8),'v':(-7,-8,9,-9),'w':(-3,-4,10,-10)}
graph_five_verts_7 = (generate_by_verts(graph_five_verts_7_by_edges),graph_five_verts_7_by_edges,[1,6,7,3])
graph_five_verts_8_by_edges = {'x':(1,2,4,5),'y':(-1,-2,3,-3),'z':(-4,6,-6,7),'v':(-7,8,-8,9),'w':(-5,-9,10,-10)}
graph_five_verts_8 = (generate_by_verts(graph_five_verts_8_by_edges),graph_five_verts_8_by_edges,[1,4,7,9])
graph_five_verts_9_by_edges = {'x':(1,-1,2,3),'y':(-2,4,-4,5),'z':(-5,6,-6,7),'v':(-7,8,-8,9),'w':(-3,-9,10,-10)}
graph_five_verts_9 = (generate_by_verts(graph_five_verts_9_by_edges),graph_five_verts_9_by_edges,[2,5,7,9])
graph_five_verts_10_by_edges = {'x':(1,-1,2,3),'y':(-2,4,-4,5),'z':(-5,6,-6,7),'v':(-7,8,9,10),'w':(-3,-8,-9,-10)}
graph_five_verts_10 = (generate_by_verts(graph_five_verts_10_by_edges),graph_five_verts_10_by_edges,[2,5,7,8])
graph_five_verts_11_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-2,-3,-6,7),'v':(-7,8,9,-9),'w':(-4,-8,10,-10)}
graph_five_verts_11 = (generate_by_verts(graph_five_verts_11_by_edges),graph_five_verts_11_by_edges,[1,6,7,8])
graph_five_verts_12_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-3,-5,-6,7),'v':(-7,8,-8,9),'w':(-4,-9,10,-10)}
graph_five_verts_12 = (generate_by_verts(graph_five_verts_12_by_edges),graph_five_verts_12_by_edges,[1,5,7,9])
graph_five_verts_13_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-5,-6,7,8),'v':(-7,-8,9,10),'w':(-3,-4,-9,-10)}
graph_five_verts_13 = (generate_by_verts(graph_five_verts_13_by_edges),graph_five_verts_13_by_edges,[1,5,7,9])
graph_five_verts_14_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-5,-6,7,8),'v':(-3,-7,9,10),'w':(-4,-8,-9,-10)}
graph_five_verts_14 = (generate_by_verts(graph_five_verts_14_by_edges),graph_five_verts_14_by_edges,[1,5,7,9])
graph_five_verts_15_by_edges = {'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-5,8,9,10),'v':(-2,-6,-8,-9),'w':(-3,-4,-7,-10)}
graph_five_verts_15 = (generate_by_verts(graph_five_verts_15_by_edges),graph_five_verts_15_by_edges,[1,3,5,8])
graph_five_verts_16_by_edges = {'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-2,-5,8,9),'v':(-3,-6,-9,10),'w':(-4,-7,-8,-10)}
graph_five_verts_16 = (generate_by_verts(graph_five_verts_16_by_edges),graph_five_verts_16_by_edges,[1,5,9,10])
graph_five_verts_17_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,-5),'z':(-3,6,7,8),'v':(-6,-7,-8,9),'w':(-4,-9,10,-10)}
graph_five_verts_17 = (generate_by_verts(graph_five_verts_17_by_edges),graph_five_verts_17_by_edges,[1,3,6,9])
graph_five_verts_18_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,-5),'z':(-3,6,7,8),'v':(-8,9,10,-10),'w':(-4,-6,-7,-9)}
graph_five_verts_18 = (generate_by_verts(graph_five_verts_18_by_edges),graph_five_verts_18_by_edges,[1,3,8,9])
graph_five_verts_19_by_edges = {'x':(1,-1,2,3),'y':(-2,4,5,6),'z':(-4,7,-7,8),'v':(-5,-8,9,-9),'w':(-3,-6,10,-10)}
graph_five_verts_19 = (generate_by_verts(graph_five_verts_19_by_edges),graph_five_verts_19_by_edges,[2,3,4,8])
graph_five_verts_20_by_edges = {'x':(1,-1,2,3),'y':(-2,4,5,6),'z':(-4,7,8,9),'v':(-5,-7,-8,-9),'w':(-3,-6,10,-10)}
graph_five_verts_20 = (generate_by_verts(graph_five_verts_20_by_edges),graph_five_verts_20_by_edges,[2,3,4,8])
graph_five_verts_21_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-2,-3,-6,7),'v':(-7,8,9,10),'w':(-4,-8,-9,-10)}
graph_five_verts_21 = (generate_by_verts(graph_five_verts_21_by_edges),graph_five_verts_21_by_edges,[1,6,7,8])
graph_five_verts_22_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-6,7,8,9),'v':(-7,-8,-9,10),'w':(-2,-3,-4,-10)}
graph_five_verts_22 = (generate_by_verts(graph_five_verts_22_by_edges),graph_five_verts_22_by_edges,[1,2,6,9])
graph_five_verts_23_by_edges = {'x':(1,2,3,4),'y':(-1,5,6,7),'z':(-5,8,9,10),'v':(-6,-8,-9,-10),'w':(-2,-3,-4,-7)}
graph_five_verts_23 = (generate_by_verts(graph_five_verts_23_by_edges),graph_five_verts_23_by_edges,[1,2,5,8])
graph_five_verts_24_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-3,-5,-6,7),'v':(-7,8,9,10),'w':(-4,-8,-9,-10)}
graph_five_verts_24 = (generate_by_verts(graph_five_verts_24_by_edges),graph_five_verts_24_by_edges,[1,4,5,7])
graph_five_verts_25_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-2,-6,7,8),'v':(-7,-8,9,10),'w':(-3,-4,-9,-10)}
graph_five_verts_25 = (generate_by_verts(graph_five_verts_25_by_edges),graph_five_verts_25_by_edges,[1,6,7,9])
graph_five_verts_26_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-2,-6,7,8),'v':(-3,-7,9,-9),'w':(-4,-8,10,-10)}
graph_five_verts_26 = (generate_by_verts(graph_five_verts_26_by_edges),graph_five_verts_26_by_edges,[1,4,6,7])
graph_five_verts_27_by_edges = {'x':(1,2,3,4),'y':(-1,-2,5,6),'z':(-5,-6,7,8),'v':(-3,-7,9,-9),'w':(-4,-8,10,-10)}
graph_five_verts_27 = (generate_by_verts(graph_five_verts_27_by_edges),graph_five_verts_27_by_edges,[1,4,5,7])
graph_five_verts_28_by_edges = {'x':(1,2,3,4),'y':(-1,5,-5,6),'z':(-6,7,8,9),'v':(-2,-7,-8,10),'w':(-3,-4,-9,-10)}
graph_five_verts_28 = (generate_by_verts(graph_five_verts_28_by_edges),graph_five_verts_28_by_edges,[1,3,6,7])

all_graphs_five_verts = [graph_five_verts_1,graph_five_verts_2,graph_five_verts_3,graph_five_verts_4,graph_five_verts_5,graph_five_verts_6,graph_five_verts_7,graph_five_verts_8,graph_five_verts_9,graph_five_verts_10,graph_five_verts_11,graph_five_verts_12,graph_five_verts_13,graph_five_verts_14,graph_five_verts_15,graph_five_verts_16,graph_five_verts_17,graph_five_verts_18,graph_five_verts_19,graph_five_verts_20,graph_five_verts_21,graph_five_verts_22,graph_five_verts_23,graph_five_verts_24,graph_five_verts_25,graph_five_verts_26,graph_five_verts_27,graph_five_verts_28]

# NEW STUFF ----------------

graph_six_verts_1 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'z', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, 10, 11, 12), 'u': (-5, -9, -10, -11), 'v': (-1, -6, -7, -12), 'w': (-2, -3, -4, -8)},[1, 6, 5, 9, 8])
graph_six_verts_2 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -10, 12, -12), 'v': (-1, -6, -7, -11), 'w': (-2, -3, -4, -8)},[1, 6, 5, 10, 8])
graph_six_verts_3 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -10, -11, 12), 'v': (-1, -6, -7, -12), 'w': (-2, -3, -4, -8)},[1, 6, 5, 10, 8])
graph_six_verts_4 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'z', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, 10, 11, 12), 'u': (-5, -9, -10, -11), 'v': (-1, -6, -7, -8), 'w': (-2, -3, -4, -12)},[1, 8, 5, 9, 12])
graph_six_verts_5 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -10, 12, -12), 'v': (-1, -6, -7, -8), 'w': (-2, -3, -4, -11)},[1, 8, 5, 10, 11])
graph_six_verts_6 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -10, -11, 12), 'v': (-1, -6, -7, -8), 'w': (-2, -3, -4, -12)},[1, 8, 5, 10, 12])
graph_six_verts_7 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'z', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, 10, 11, 12), 'u': (-5, -6, -9, -10), 'v': (-1, -7, -11, -12), 'w': (-2, -3, -4, -8)},[1, 7, 5, 9, 8])
graph_six_verts_8 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'z', -11: 'v', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, 12, -12), 'v': (-1, -7, -10, -11), 'w': (-2, -3, -4, -8)},[1, 7, 5, 8, 10])
graph_six_verts_9 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, 12), 'v': (-1, -7, -11, -12), 'w': (-2, -3, -4, -8)},[1, 7, 5, 10, 8])
graph_six_verts_10 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, -11), 'v': (-1, -7, 12, -12), 'w': (-2, -3, -4, -8)},[1, 7, 5, 10, 8])
graph_six_verts_11 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, 12), 'v': (-1, -7, -8, -12), 'w': (-2, -3, -4, -11)},[1, 8, 5, 10, 11])
graph_six_verts_12 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, 12), 'v': (-1, -7, -8, -11), 'w': (-2, -3, -4, -12)},[1, 8, 5, 10, 12])
graph_six_verts_13 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, -11), 'v': (-1, -7, -8, 12), 'w': (-2, -3, -4, -12)},[1, 8, 5, 10, 12])
graph_six_verts_14 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -7, 12), 'v': (-1, -8, -10, -12), 'w': (-2, -3, -4, -11)},[1, 8, 5, 10, 11])
graph_six_verts_15 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -7, -10), 'v': (-1, -8, 12, -12), 'w': (-2, -3, -4, -11)},[1, 8, 5, 10, 11])
graph_six_verts_16 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -9, -10, 12), 'v': (-1, -7, -11, -12), 'w': (-2, -3, -4, -8)},[1, 7, 5, 9, 8])
graph_six_verts_17 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, 11, -11, 12), 'v': (-1, -7, -10, -12), 'w': (-2, -3, -4, -8)},[1, 7, 5, 6, 8])
graph_six_verts_18 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -10, 11, -11), 'v': (-1, -7, 12, -12), 'w': (-2, -3, -4, -8)},[1, 7, 5, 10, 8])
graph_six_verts_19 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -9, -10, 12), 'v': (-1, -7, -8, -12), 'w': (-2, -3, -4, -11)},[1, 8, 5, 9, 11])
graph_six_verts_20 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -9, -10, -11), 'v': (-1, -7, -8, 12), 'w': (-2, -3, -4, -12)},[1, 8, 5, 9, 12])
graph_six_verts_21 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, 11, -11, 12), 'v': (-1, -7, -8, -12), 'w': (-2, -3, -4, -10)},[1, 8, 5, 10, 6])
graph_six_verts_22 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -10, 11, -11), 'v': (-1, -7, -8, 12), 'w': (-2, -3, -4, -12)},[1, 8, 5, 10, 12])
graph_six_verts_23 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, 11, -11), 'v': (-1, -10, 12, -12), 'w': (-2, -3, -4, -8)},[1, 10, 5, 6, 8])
graph_six_verts_24 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -7, -9, -10), 'v': (-1, -8, -11, 12), 'w': (-2, -3, -4, -12)},[1, 8, 5, 9, 12])
graph_six_verts_25 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, 11, -11), 'v': (-1, -8, 12, -12), 'w': (-2, -3, -4, -10)},[1, 8, 5, 10, 6])
graph_six_verts_26 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, 11, -11), 'v': (-1, -8, -10, 12), 'w': (-2, -3, -4, -12)},[1, 8, 5, 6, 12])
graph_six_verts_27 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'v', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, -10, 11), 'v': (-1, -8, -11, 12), 'w': (-2, -3, -4, -12)},[1, 8, 5, 10, 12])
graph_six_verts_28 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-6, -9, 11, -11), 'v': (-1, -10, 12, -12), 'w': (-2, -3, -4, -7)},[1, 10, 9, 6, 7])
graph_six_verts_29 = ({1: 'x', -1: 'v', 2: 'x', -2: 'w', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-6, -9, 11, -11), 'v': (-1, -7, -10, 12), 'w': (-2, -3, -4, -12)},[1, 7, 6, 9, 12])
graph_six_verts_30 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -10, 12, -12), 'v': (-1, -2, -6, -11), 'w': (-3, -4, -7, -8)},[1, 6, 5, 10, 8])
graph_six_verts_31 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -10, -11, 12), 'v': (-1, -2, -6, -12), 'w': (-3, -4, -7, -8)},[1, 6, 5, 10, 8])
graph_six_verts_32 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'z', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, 10, 11, 12), 'u': (-5, -6, -9, -10), 'v': (-1, -2, -11, -12), 'w': (-3, -4, -7, -8)},[1, 11, 9, 5, 8])
graph_six_verts_33 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'z', -11: 'v', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, 12, -12), 'v': (-1, -2, -10, -11), 'w': (-3, -4, -7, -8)},[1, 10, 3, 8, 5])
graph_six_verts_34 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, 12), 'v': (-1, -2, -11, -12), 'w': (-3, -4, -7, -8)},[1, 11, 10, 5, 8])
graph_six_verts_35 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'z', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, 10, 11, 12), 'u': (-5, -6, -9, -10), 'v': (-1, -2, -7, -11), 'w': (-3, -4, -8, -12)},[1, 7, 5, 9, 12])
graph_six_verts_36 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'z', -11: 'w', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, 12, -12), 'v': (-1, -2, -7, -10), 'w': (-3, -4, -8, -11)},[1, 7, 5, 8, 11])
graph_six_verts_37 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, 12), 'v': (-1, -2, -7, -12), 'w': (-3, -4, -8, -11)},[1, 7, 5, 10, 11])
graph_six_verts_38 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-5, -6, -10, -11), 'v': (-1, -2, -7, 12), 'w': (-3, -4, -8, -12)},[1, 7, 5, 10, 8])
graph_six_verts_39 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -9, -10, 12), 'v': (-1, -2, -11, -12), 'w': (-3, -4, -7, -8)},[1, 11, 5, 6, 8])
graph_six_verts_40 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, 11, -11, 12), 'v': (-1, -2, -10, -12), 'w': (-3, -4, -7, -8)},[1, 10, 5, 6, 8])
graph_six_verts_41 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -10, 11, -11), 'v': (-1, -2, 12, -12), 'w': (-3, -4, -7, -8)},[1, 3, 8, 5, 10])
graph_six_verts_42 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'z', -11: 'w', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -9, 12, -12), 'v': (-1, -2, -7, -10), 'w': (-3, -4, -8, -11)},[1, 7, 5, 9, 11])
graph_six_verts_43 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -9, -10, 12), 'v': (-1, -2, -7, -12), 'w': (-3, -4, -8, -11)},[1, 7, 5, 9, 11])
graph_six_verts_44 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, 11, -11, 12), 'v': (-1, -2, -7, -12), 'w': (-3, -4, -8, -10)},[1, 7, 5, 10, 6])
graph_six_verts_45 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -10, 11, -11), 'v': (-1, -2, -7, 12), 'w': (-3, -4, -8, -12)},[1, 7, 5, 10, 8])
graph_six_verts_46 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-6, -7, -9, -10), 'v': (-1, -2, -11, 12), 'w': (-3, -4, -8, -12)},[1, 11, 5, 6, 8])
graph_six_verts_47 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, 11, -11), 'v': (-1, -2, 12, -12), 'w': (-3, -4, -8, -10)},[1, 3, 8, 5, 6])
graph_six_verts_48 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, 11, -11), 'v': (-1, -2, -10, 12), 'w': (-3, -4, -8, -12)},[1, 10, 5, 6, 8])
graph_six_verts_49 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'w', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, -10, 11), 'v': (-1, -2, 12, -12), 'w': (-3, -4, -8, -11)},[1, 3, 8, 5, 10])
graph_six_verts_50 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'u', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'v', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-6, -7, -10, 11), 'v': (-1, -2, -11, 12), 'w': (-3, -4, -8, -12)},[1, 11, 6, 5, 8])
graph_six_verts_51 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-9, 11, -11, 12), 'v': (-1, -2, -6, -12), 'w': (-3, -4, -7, -10)},[1, 6, 7, 10, 9])
graph_six_verts_52 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-6, -9, 11, -11), 'v': (-1, -2, 12, -12), 'w': (-3, -4, -7, -10)},[1, 3, 7, 6, 9])
graph_six_verts_53 = ({1: 'x', -1: 'v', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'u', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-6, -9, 11, -11), 'v': (-1, -2, -10, 12), 'w': (-3, -4, -7, -12)},[1, 10, 9, 6, 7])
graph_six_verts_54 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'z', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, 10, 11, 12), 'u': (-1, -5, -9, -10), 'v': (-2, -6, -7, -11), 'w': (-3, -4, -8, -12)},[1, 5, 6, 11, 12])
graph_six_verts_55 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'z', -11: 'w', 12: 'u', -12: 'u'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-1, -5, 12, -12), 'v': (-2, -6, -7, -10), 'w': (-3, -4, -8, -11)},[1, 5, 6, 10, 11])
graph_six_verts_56 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-1, -5, -10, 12), 'v': (-2, -6, -7, -12), 'w': (-3, -4, -8, -11)},[1, 5, 6, 8, 11])
graph_six_verts_57 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'u', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'z', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (9, -9, 10, 11), 'u': (-1, -5, -10, -11), 'v': (-2, -6, -7, 12), 'w': (-3, -4, -8, -12)},[1, 5, 6, 12, 10])
graph_six_verts_58 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-1, -9, -10, 12), 'v': (-2, -6, -7, -12), 'w': (-3, -4, -8, -11)},[1, 9, 5, 6, 8])
graph_six_verts_59 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'z', -11: 'v', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-1, -9, -10, 12), 'v': (-2, -6, -7, -11), 'w': (-3, -4, -8, -12)},[1, 9, 5, 6, 8])
graph_six_verts_60 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, 11, -11, 12), 'v': (-2, -6, -7, -12), 'w': (-3, -4, -8, -10)},[1, 12, 6, 5, 10])
graph_six_verts_61 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, 11, -11, 12), 'v': (-2, -6, -7, -10), 'w': (-3, -4, -8, -12)},[1, 12, 8, 5, 10])
graph_six_verts_62 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -10, 11, -11), 'v': (-2, -6, -7, 12), 'w': (-3, -4, -8, -12)},[1, 10, 5, 6, 12])
graph_six_verts_63 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'z', -11: 'w', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-1, -6, -9, 12), 'v': (-2, -7, -10, -12), 'w': (-3, -4, -8, -11)},[1, 6, 5, 10, 11])
graph_six_verts_64 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'z', -11: 'v', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, 10, 11), 'u': (-1, -6, -9, 12), 'v': (-2, -7, -10, -11), 'w': (-3, -4, -8, -12)},[1, 6, 5, 10, 8])
graph_six_verts_65 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'v', 12: 'u', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, 12), 'v': (-2, -7, -11, -12), 'w': (-3, -4, -8, -10)},[1, 6, 5, 10, 7])
graph_six_verts_66 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, -11), 'v': (-2, -7, 12, -12), 'w': (-3, -4, -8, -10)},[1, 6, 5, 10, 7])
graph_six_verts_67 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'v', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, 12), 'v': (-2, -7, -10, -11), 'w': (-3, -4, -8, -12)},[1, 6, 5, 10, 8])
graph_six_verts_68 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, -11), 'v': (-2, -7, -10, 12), 'w': (-3, -4, -8, -12)},[1, 6, 5, 10, 12])
graph_six_verts_69 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'v', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, 12), 'v': (-2, -7, -8, -11), 'w': (-3, -4, -10, -12)},[1, 6, 5, 10, 8])
graph_six_verts_70 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, -11), 'v': (-2, -7, -8, 12), 'w': (-3, -4, -10, -12)},[1, 6, 5, 10, 12])
graph_six_verts_71 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, 11, -11), 'v': (-2, -7, -8, -10), 'w': (-3, -4, 12, -12)},[1, 6, 5, 10, 3])
graph_six_verts_72 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'z', -10: 'u', 11: 'u', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, 9, -9, 10), 'u': (-1, -6, -10, 11), 'v': (-2, -7, -8, -11), 'w': (-3, -4, 12, -12)},[1, 6, 5, 8, 3])
graph_six_verts_73 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'z', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'w', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, -6, 9, 10), 'u': (-1, -7, -9, 11), 'v': (-2, -8, -10, 12), 'w': (-3, -4, -11, -12)},[1, 7, 5, 10, 12])
graph_six_verts_74 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'z', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, -6, 9, 10), 'u': (-1, -7, -9, 11), 'v': (-2, -8, -10, -11), 'w': (-3, -4, 12, -12)},[1, 7, 5, 10, 3])
graph_six_verts_75 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'z', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'u', -10: 'v', 11: 'u', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, -6, 9, -9), 'u': (-1, -7, 10, 11), 'v': (-2, -8, -10, -11), 'w': (-3, -4, 12, -12)},[1, 7, 5, 8, 3])
graph_six_verts_76 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'z', 7: 'y', -7: 'u', 8: 'y', -8: 'v', 9: 'z', -9: 'z', 10: 'u', -10: 'u', 11: 'v', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-5, -6, 9, -9), 'u': (-1, -7, 10, -10), 'v': (-2, -8, 11, -11), 'w': (-3, -4, 12, -12)},[1, 7, 5, 8, 3])
graph_six_verts_77 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -9, 11, -11), 'v': (-2, -6, 12, -12), 'w': (-3, -4, -7, -10)},[1, 9, 10, 7, 6])
graph_six_verts_78 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -9, 11, -11), 'v': (-2, -6, -10, 12), 'w': (-3, -4, -7, -12)},[1, 9, 10, 6, 7])
graph_six_verts_79 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'u', -11: 'w', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -9, -10, 11), 'v': (-2, -6, 12, -12), 'w': (-3, -4, -7, -11)},[1, 9, 11, 7, 6])
graph_six_verts_80 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'u', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -9, 11, -11), 'v': (-2, -6, -7, -10), 'w': (-3, -4, 12, -12)},[1, 9, 10, 6, 3])
graph_six_verts_81 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'u', 11: 'u', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -9, -10, 11), 'v': (-2, -6, -7, -11), 'w': (-3, -4, 12, -12)},[1, 9, 11, 6, 3])
graph_six_verts_82 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'w', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -6, -9, 11), 'v': (-2, -7, -10, 12), 'w': (-3, -4, -11, -12)},[1, 6, 7, 10, 12])
graph_six_verts_83 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'z', -10: 'v', 11: 'u', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (8, -8, 9, 10), 'u': (-1, -6, -9, 11), 'v': (-2, -7, -10, -11), 'w': (-3, -4, 12, -12)},[1, 6, 7, 10, 3])
graph_six_verts_84 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'z', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'w', 10: 'u', -10: 'u', 11: 'u', -11: 'v', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-6, 8, -8, 9), 'u': (-1, 10, -10, 11), 'v': (-2, -11, 12, -12), 'w': (-3, -4, -7, -9)},[1, 11, 3, 7, 6])
graph_six_verts_85 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'z', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'v', 10: 'u', -10: 'u', 11: 'u', -11: 'w', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-6, 8, -8, 9), 'u': (-1, 10, -10, 11), 'v': (-2, -9, 12, -12), 'w': (-3, -4, -7, -11)},[1, 11, 7, 6, 9])
graph_six_verts_86 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'z', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'v', 10: 'u', -10: 'u', 11: 'u', -11: 'w', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-6, 8, -8, 9), 'u': (-1, 10, -10, 11), 'v': (-2, -7, -9, 12), 'w': (-3, -4, -11, -12)},[1, 11, 12, 7, 6])
graph_six_verts_87 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'z', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'v', 10: 'u', -10: 'u', 11: 'u', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-6, 8, -8, 9), 'u': (-1, 10, -10, 11), 'v': (-2, -7, -9, -11), 'w': (-3, -4, 12, -12)},[1, 11, 7, 6, 3])
graph_six_verts_88 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'z', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'u', -10: 'v', 11: 'u', -11: 'w', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-6, 8, -8, 9), 'u': (-1, -9, 10, 11), 'v': (-2, -7, -10, 12), 'w': (-3, -4, -11, -12)},[1, 9, 6, 7, 12])
graph_six_verts_89 = ({1: 'x', -1: 'u', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'z', 7: 'y', -7: 'v', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'u', -10: 'u', 11: 'v', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-6, 8, -8, 9), 'u': (-1, -9, 10, -10), 'v': (-2, -7, 11, -11), 'w': (-3, -4, 12, -12)},[1, 9, 6, 7, 3])
graph_six_verts_90 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'v', 10: 'z', -10: 'w', 11: 'u', -11: 'v', 12: 'u', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-1, -5, 9, 10), 'u': (-2, -6, 11, 12), 'v': (-3, -7, -9, -11), 'w': (-4, -8, -10, -12)},[1, 5, 6, 11, 12])
graph_six_verts_91 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'v', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-1, -5, 9, 10), 'u': (-2, -6, 11, -11), 'v': (-3, -7, -9, 12), 'w': (-4, -8, -10, -12)},[1, 5, 6, 7, 12])
graph_six_verts_92 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'z', 6: 'y', -6: 'u', 7: 'y', -7: 'v', 8: 'y', -8: 'w', 9: 'z', -9: 'z', 10: 'u', -10: 'u', 11: 'v', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, 6, 7, 8), 'z': (-1, -5, 9, -9), 'u': (-2, -6, 10, -10), 'v': (-3, -7, 11, -11), 'w': (-4, -8, 12, -12)},[1, 5, 6, 7, 8])
graph_six_verts_93 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'u', 9: 'z', -9: 'v', 10: 'z', -10: 'w', 11: 'u', -11: 'u', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-1, 8, 9, 10), 'u': (-2, -8, 11, -11), 'v': (-3, -6, -9, 12), 'w': (-4, -7, -10, -12)},[1, 8, 9, 6, 7])
graph_six_verts_94 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'w', 10: 'u', -10: 'u', 11: 'u', -11: 'w', 12: 'v', -12: 'v'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-1, 8, -8, 9), 'u': (-2, 10, -10, 11), 'v': (-3, -6, 12, -12), 'w': (-4, -7, -9, -11)},[1, 9, 7, 6, 11])
graph_six_verts_95 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'w', 10: 'u', -10: 'u', 11: 'u', -11: 'v', 12: 'v', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-1, 8, -8, 9), 'u': (-2, 10, -10, 11), 'v': (-3, -6, -11, 12), 'w': (-4, -7, -9, -12)},[1, 9, 7, 6, 11])
graph_six_verts_96 = ({1: 'x', -1: 'z', 2: 'x', -2: 'u', 3: 'x', -3: 'v', 4: 'x', -4: 'w', 5: 'y', -5: 'y', 6: 'y', -6: 'v', 7: 'y', -7: 'w', 8: 'z', -8: 'z', 9: 'z', -9: 'u', 10: 'u', -10: 'u', 11: 'v', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, 2, 3, 4), 'y': (5, -5, 6, 7), 'z': (-1, 8, -8, 9), 'u': (-2, -9, 10, -10), 'v': (-3, -6, 11, -11), 'w': (-4, -7, 12, -12)},[1, 9, 3, 6, 7])
graph_six_verts_97 = ({1: 'x', -1: 'x', 2: 'x', -2: 'v', 3: 'x', -3: 'w', 4: 'y', -4: 'y', 5: 'y', -5: 'u', 6: 'y', -6: 'w', 7: 'z', -7: 'z', 8: 'z', -8: 'u', 9: 'z', -9: 'v', 10: 'u', -10: 'u', 11: 'v', -11: 'v', 12: 'w', -12: 'w'},{'x': (1, -1, 2, 3), 'y': (4, -4, 5, 6), 'z': (7, -7, 8, 9), 'u': (-5, -8, 10, -10), 'v': (-2, -9, 11, -11), 'w': (-3, -6, 12, -12)},[2, 9, 8, 5, 6])

all_graphs_six_verts = [graph_six_verts_1, graph_six_verts_2, graph_six_verts_3, graph_six_verts_4, graph_six_verts_5, graph_six_verts_6, graph_six_verts_7, graph_six_verts_8, graph_six_verts_9, graph_six_verts_10, graph_six_verts_11, graph_six_verts_12, graph_six_verts_13, graph_six_verts_14, graph_six_verts_15, graph_six_verts_16, graph_six_verts_17, graph_six_verts_18, graph_six_verts_19, graph_six_verts_20, graph_six_verts_21, graph_six_verts_22, graph_six_verts_23, graph_six_verts_24, graph_six_verts_25, graph_six_verts_26, graph_six_verts_27, graph_six_verts_28, graph_six_verts_29, graph_six_verts_30, graph_six_verts_31, graph_six_verts_32, graph_six_verts_33, graph_six_verts_34, graph_six_verts_35, graph_six_verts_36, graph_six_verts_37, graph_six_verts_38, graph_six_verts_39, graph_six_verts_40, graph_six_verts_41, graph_six_verts_42, graph_six_verts_43, graph_six_verts_44, graph_six_verts_45, graph_six_verts_46, graph_six_verts_47, graph_six_verts_48, graph_six_verts_49, graph_six_verts_50, graph_six_verts_51, graph_six_verts_52, graph_six_verts_53, graph_six_verts_54, graph_six_verts_55, graph_six_verts_56, graph_six_verts_57, graph_six_verts_58, graph_six_verts_59, graph_six_verts_60, graph_six_verts_61, graph_six_verts_62, graph_six_verts_63, graph_six_verts_64, graph_six_verts_65, graph_six_verts_66, graph_six_verts_67, graph_six_verts_68, graph_six_verts_69, graph_six_verts_70, graph_six_verts_71, graph_six_verts_72, graph_six_verts_73, graph_six_verts_74, graph_six_verts_75, graph_six_verts_76, graph_six_verts_77, graph_six_verts_78, graph_six_verts_79, graph_six_verts_80, graph_six_verts_81, graph_six_verts_82, graph_six_verts_83, graph_six_verts_84, graph_six_verts_85, graph_six_verts_86, graph_six_verts_87, graph_six_verts_88, graph_six_verts_89, graph_six_verts_90, graph_six_verts_91, graph_six_verts_92, graph_six_verts_93, graph_six_verts_94, graph_six_verts_95, graph_six_verts_96, graph_six_verts_97]
# --------------

arg_list = sys.argv
graph_num = int(arg_list[1])
print("running graph number: " + str(graph_num))
curr_graph = all_graphs_six_verts[graph_num-1]

all_disks = create_disks(curr_graph) #dfs part 
unique_disks = remove_duplicate_disks(all_disks)
acyclic_fake_surfaces, all_partitions, acyclic_partitions = find_acyclic_fake_surfaces(unique_disks, curr_graph)
unique_acyclic_fake_surfaces = remove_duplicate_fake_surfaces(acyclic_fake_surfaces)

# print(unique_acyclic_fake_surfaces)
# print("Number of unique acyclic fake surfaces:",sum([len(unique_acyclic_fake_surfaces[key]) for key in unique_acyclic_fake_surfaces.keys()]))
# print("Unique acyclic fake surface partitions:",{key:len(unique_acyclic_fake_surfaces[key]) for key in unique_acyclic_fake_surfaces.keys()})
with open('/home/lucasfagan/fakesurfaces_complexity6/unique_surfaces_graph_six_verts_'+str(graph_num)+'.txt', 'w') as f:
    f.write("Number of unique acyclic fake surfaces: "+str(sum([len(unique_acyclic_fake_surfaces[key]) for key in unique_acyclic_fake_surfaces.keys()]))+'\n')
    f.write(str(unique_acyclic_fake_surfaces)+'\n')

fake_surfaces_no_dict = []
for key in unique_acyclic_fake_surfaces.keys():
    fake_surfaces_no_dict.extend(unique_acyclic_fake_surfaces[key])
write_to_file_for_sage(fake_surfaces_no_dict, curr_graph, graph_num)

