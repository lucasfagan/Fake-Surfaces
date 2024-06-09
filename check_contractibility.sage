# sample data from complexity 2 to demonstrate formatting
acyclic_surface = [[4,2,-1,-1,-2,3,2,1,-3],[2,-3],[4]]
maximal_tree = [2]

def is_contractible(acyclic_surface, maximal_tree):
    # remove the elements in the maximal tree
    collapsed_surface = [[edge for edge in disk if (edge not in maximal_tree and -1*edge not in maximal_tree)] for disk in acyclic_surface]
    # for conveinence, reduce the remaining edges to be labeled with 1 through n+1
    for removed_edge in sorted(maximal_tree, reverse=True):
        collapsed_surface = [[edge-1 if edge>removed_edge else edge for edge in disk] for disk in collapsed_surface]
        collapsed_surface = [[edge+1 if edge<-1*removed_edge else edge for edge in disk] for disk in collapsed_surface] 

    # create a free group with one generator for each remaining edge
    F = FreeGroup(len(maximal_tree)+2,'a')

    group_relations = []
    
    # convert the disks to relations 
    for disk in collapsed_surface:
        relation = F.one()  # Initialize with identity
        for edge in disk:
            if edge < 0:
                relation *= (F.gen(-edge-1))**-1  # Inverse of generator
            else:
                relation *= F.gen(edge-1)  # Multiply with generator

        # Append the relation
        group_relations.append(relation)

    # Create quotient group
    G = F.quotient(group_relations) 
    # Simplify the group
    G = G.simplified()

    # Check if group is trivial
    if G.order()!=1:
        return False
    return True
