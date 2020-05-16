def get_lt_adj_matrix_entry(i, j, adj_matrix):
    # lt = left-triangular
    if len(adj_matrix[i]) <= j:
        return adj_matrix[j][i]

    return adj_matrix[i][j]


def create_adj_matrix_for_comms(adj_matrix, partition, starting_node):
    comm_adj_matrix = []
    buffer_edge_vecs = [[0.0] for _ in range(starting_node)]
    for i, community in enumerate(partition):
        comm_edge_vecs = []
        for j, neighbor_comm in enumerate(partition):
            if j > i:
                break

            edge_vecs = []
            for u, v in product(community, neighbor_comm):
                edge_vec = get_lt_adj_matrix_entry(u, v, adj_matrix)
                edge_vecs.append(edge_vec)

            avg_edge_vec = list(map(mean, zip(*edge_vecs)))
            comm_edge_vecs.append(avg_edge_vec)

        comm_adj_matrix.append(buffer_edge_vecs + comm_edge_vecs) # QUESTION: What are these buffer vecs for?

    return comm_adj_matrix


def create_comm_graphs(partition, adj_matrix):
    for community in partition:
        comm_adj_matrix = []
        for i, src in enumerate(community):
            edge_vecs = []
            for j, targ in enumerate(community):
                if j > i:
                    break

                edge_vec = get_lt_adj_matrix_entry(src, targ, adj_matrix)
                edge_vecs.append(edge_vec)

            comm_adj_matrix.append(edge_vecs)

        yield comm_adj_matrix
