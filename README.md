# network_robustness_max_vulnerability
Calculation of network Robustness and maximum vulnerability in R

**Network Robustness**

Robustness, defined as the capacity of a network to maintain connectivity under node removal, was assessed using two strategies: random removal and targeted removal with random tie-breaking. In random removal, nodes were selected randomly, while targeted removal sequentially removed nodes based on degree centrality, targeting the highest-degree nodes first. 
Ties (nodes with equal degrees) were resolved by randomly selecting one node, introducing variability into the targeted strategy. See: AUC_random_targeted_node_removal.R

Group-specific robustness was also analyzed for predefined microbial groups (e.g., fungi and protists). Nodes corresponding to these groups were sequentially removed, with the removal order randomized across 100 iterations. 
The largest connected component size was recalculated after each removal, and AUC metrics were derived for each group, enabling assessment of their contributions to overall network connectivity. See: AUC_targeted_trophic_group_removal.R

**Network Vulnerability Analysis**

Network vulnerability, reflecting the relative reduction in global efficiency after node removal, was evaluated for each network and microbial group. Briefly, global network efficiency (E) was calculated for each network and recalculated after node removal. 
The removal process involved deleting the node and all its associated edges using the igraph function delete_vertices(). 
Vulnerabilities were computed for each node within a group, providing a detailed assessment of how individual nodes contribute to the network's global efficiency. 
Maximum vulnerability, representing the most critical node, and the mean and standard deviation of vulnerabilities across all nodes were computed for each network. See: max_vulnerability.R
