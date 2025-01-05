#The following script was generated to calculate network robustnes (AUC) after randomized removal of specific network groups, plotting and statistical testing
#Here, we removed higher trophic groups (fungi and protists) separately from the network to calculate its robustness
#networks are imported as tab generated files

#set working directory
setwd("/Users/networks") #this corresponds to the location of your network's files

#load package
library(igraph)

#import network edge data
edge_lists.T1 <- as.data.frame(read.table("allT1nonTrare002_0.96_edge_attribute_idirect.txt", header = T, as.is = T, sep="\t"))
edge_lists.T3 <- as.data.frame(read.table("allT3nonTrare002_0.96_edge_attribute_idirect.txt", header = T, as.is = T, sep="\t"))
edge_lists.T4 <- as.data.frame(read.table("allT4nonTrare002_0.94_edge_attribute_idirect.txt", header = T, as.is = T, sep="\t"))

#declare networks
network1 <- graph_from_edgelist(as.matrix(edge_lists.T1[, c("taxa1", "taxa2")]), directed = FALSE)
network2 <-graph_from_edgelist(as.matrix(edge_lists.T3[, c("taxa1", "taxa2")]), directed = FALSE)
network3 <-graph_from_edgelist(as.matrix(edge_lists.T4[, c("taxa1", "taxa2")]), directed = FALSE)


#attempt specific removal by node list or kingdom


comm_results.T1 <- as.data.frame(read.table("allT1_nonTrare002_0.96_node_attribute_idirect.txt", header = T, as.is = T, sep="\t"))
comm_results.T3 <- as.data.frame(read.table("allT3nonTrare002_0.96_node_attribute_idirect.txt", header = T, as.is = T, sep="\t"))
comm_results.T4 <- as.data.frame(read.table("allT4nonTrare002_0.94_node_attribute_idirect.txt", header = T, as.is = T, sep="\t"))


# Add kingdom information to graph nodes
add_kingdom_to_network <- function(network, comm_results) {
  V(network)$kingdom <- comm_results$kingdom[match(V(network)$name, comm_results$taxon)]
  return(network)
}

# Add kingdom information to all networks
network1 <- add_kingdom_to_network(network1, comm_results.T1)
network2 <- add_kingdom_to_network(network2, comm_results.T3)
network3 <- add_kingdom_to_network(network3, comm_results.T4)

# Targeted removal function
targeted_removal_general <- function(graph, target_nodes) {
  graph_copy <- graph  # Copy the graph
  component_sizes <- numeric(vcount(graph_copy))  # Initialize size tracker
  
  for (node in target_nodes) {
    if (node %in% V(graph_copy)$name) {  # Check if node exists
      graph_copy <- delete_vertices(graph_copy, node)  # Remove node
    }
    if (vcount(graph_copy) > 0) {
      comp <- components(graph_copy)
      component_sizes[vcount(graph_copy)] <- max(comp$csize)  # Largest component size
    } else {
      component_sizes[vcount(graph_copy)] <- 0  # No nodes left
    }
  }
  return(component_sizes)
}

# AUC calculation function with shuffled removal order
calculate_auc <- function(graph, target_nodes, iterations = 100) {
  if (length(target_nodes) == 0) {  # Handle cases with no target nodes
    warning("No target nodes available for removal.")
    return(rep(NA, iterations))
  }
  
  auc_values <- numeric(iterations)
  
  for (i in seq_len(iterations)) {
    shuffled_targets <- sample(target_nodes)  # Shuffle target node order
    result <- targeted_removal_general(graph, shuffled_targets)  # Perform removal
    auc_values[i] <- (sum(result) / length(result)) / ((vcount(graph) + 1) / 2)  # Calculate normalized AUC
  }
  return(auc_values)
}


# Define groups and networks
groups <- c("Fungi", "Protists")
networks <- list(Network1 = network1, Network2 = network2, Network3 = network3)

# Initialize an empty data frame to store raw AUC values
results_long <- data.frame(
  Network = character(),
  Removed_Group = character(),
  AUC = numeric()
)

# Iterate over networks and groups
for (network_name in names(networks)) {
  for (group in groups) {
    if (group %in% unique(V(networks[[network_name]])$kingdom)) {  # Check if group exists
      auc_values <- calculate_auc(networks[[network_name]], V(networks[[network_name]])$name[V(networks[[network_name]])$kingdom == group])
      results_long <- rbind(
        results_long,
        data.frame(Network = rep(network_name, length(auc_values)),
                   Removed_Group = rep(group, length(auc_values)),
                   AUC = auc_values)
      )
    }
  }
}

# View the raw data frame
print(results_long)

library(dplyr)
# Load ggplot2
library(ggplot2)

# Aggregate data to calculate mean and SD for each group and network
results_summary <- results_long %>%
  group_by(Network, Removed_Group) %>%
  summarize(
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    .groups = "drop"
  )
results_summary

write.csv(results_summary, "removal_group_121024.csv")

removed <- ggplot(results_summary, aes(x = Removed_Group, y = Mean_AUC, fill = Network)) +
  # Bar plot
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  # Overlay raw data points with black border
  geom_jitter(
    data = results_long, 
    aes(x = Removed_Group, y = AUC, fill = Network),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
    alpha = 0.8, size = 3, shape = 21, stroke = 0.5, color = "black"
  ) +
  # Error bars
  geom_errorbar(
    aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC),
    position = position_dodge(width = 0.9), width = 0.2, size = 0.5 
  ) +
  # Customize the appearance
  labs(
    title = "AUC Across Groups and Networks",
    x = "Removed Group",
    y = "AUC",
    fill = "Network",
    color = "Network"
  ) +
  scale_y_continuous(limits = c(0, 0.4))+
  scale_fill_manual(
    values = c("Network1" = "#009E73", "Network2" = "#CC79A7", "Network3" = "#D55E00")) +
  scale_color_manual(
    values = c("Network1" = "#009E73", "Network2" = "#CC79A7", "Network3" = "#D55E00"))+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

removed
ggsave(removed, filename = "removed_121024.pdf", height = 4, width = 5)

#check normality
normality_results <- results_long %>%
  group_by(Network, Removed_Group) %>%
  summarise(
    p_value = shapiro.test(AUC)$p.value
  ) %>%
  mutate(Normality = ifelse(p_value > 0.05, "Normal", "Non-Normal"))

print(normality_results)

# Perform ANOVA
anova_result <- aov(AUC ~ Removed_Group * Network, data = results_long)
summary(anova_result)

# Post-hoc Tukey test
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(AUC ~ interaction(Removed_Group, Network), data = results_long)
print(kruskal_result)

# Post-hoc Dunn's test
library(FSA)
dunn_result <- dunnTest(AUC ~ interaction(Removed_Group, Network), data = results_long, method = "fdr")
print(dunn_result)




