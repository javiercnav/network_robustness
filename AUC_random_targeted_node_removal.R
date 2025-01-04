#The following script was generated to calculate network robustnes (AUC) under random and targeted node removal, plotting and statistical testing
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



#load functions
#single random node removal , required for multiple randomn removal with iterations
random_removal <- function(graph) {
  graph <- simplify(graph)                        # Ensure no duplicate edges or loops
  component_sizes <- numeric(vcount(graph))       # Initialize vector for component sizes
  
  for (i in seq_len(vcount(graph))) {
    nodes <- V(graph)                             # Get current nodes dynamically
    node_to_remove <- sample(nodes, 1)           # Randomly select one node to remove
    graph <- delete_vertices(graph, node_to_remove) # Remove the selected node
    if (vcount(graph) > 0) {                      # Check if nodes remain
      comp <- components(graph)
      component_sizes[i] <- max(comp$csize)       # Size of the largest connected component
    } else {
      component_sizes[i] <- 0                    # No nodes left, component size is 0
    }
  }
  return(component_sizes)
}

# Function to perform multiple random removal iterations and normalization by network size
random_removal_iterations <- function(graph, iterations = 1000) {
  auc_values <- numeric(iterations)  # Store AUC for each iteration
  
  for (i in seq_len(iterations)) {
    component_sizes <- random_removal(graph)  # Perform random node removal
    auc_values[i] <- (sum(component_sizes) / length(component_sizes)) / ((vcount(graph) + 1) / 2)  # Calculate normalized AUC
  }
  
  return(auc_values)  # Return AUCs for all iterations
}


# Targeted Node Removal required for iterative function
targeted_removal <- function(graph) {
  graph <- simplify(graph)                        # Ensure no duplicate edges or loops
  component_sizes <- numeric(vcount(graph))       # Initialize vector for component sizes
  
  for (i in seq_along(component_sizes)) {
    if (vcount(graph) > 0) {
      degree_list <- degree(graph)                # Calculate node degrees
      node_to_remove <- which.max(degree_list)    # Select the highest-degree node
      graph <- delete_vertices(graph, node_to_remove) # Remove the node
      if (vcount(graph) > 0) {
        comp <- components(graph)
        component_sizes[i] <- max(comp$csize)     # Size of the largest connected component
      } else {
        component_sizes[i] <- 0                  # No nodes left, component size is 0
      }
    }
  }
  return(component_sizes)
}


# Targeted removal function with random tie-breaking (certain nodes may have same degree and tie breaking is needed)
targeted_removal_random_ties <- function(graph, iterations = 1000) {
  auc_values <- numeric(iterations)
  
  for (i in seq_len(iterations)) {
    graph_copy <- graph
    component_sizes <- numeric(vcount(graph_copy))
    
    while (vcount(graph_copy) > 0) {
      degree_list <- degree(graph_copy)
      max_degree_nodes <- which(degree_list == max(degree_list))  # Nodes with max degree
      node_to_remove <- sample(max_degree_nodes, 1)  # Randomly select one if ties exist
      graph_copy <- delete_vertices(graph_copy, node_to_remove)
      if (vcount(graph_copy) > 0) {
        comp <- components(graph_copy)
        component_sizes[vcount(graph_copy)] <- max(comp$csize)
      } else {
        component_sizes[vcount(graph_copy)] <- 0
      }
    }
    auc_values[i] <- (sum(component_sizes) / length(component_sizes)) / ((vcount(graph) + 1) / 2)  # Calculate normalized AUC
  }
  
  return(auc_values)
}

#perform calculations

# Run with increased iterations
random_auc_network1 <- random_removal_iterations(network1, iterations = 1000)
random_auc_network2 <- random_removal_iterations(network2, iterations = 1000)
random_auc_network3 <- random_removal_iterations(network3, iterations = 1000)

# Run targeted removal with random tie-breaking
targeted_auc_network1 <- targeted_removal_random_ties(network1, iterations = 1000)
targeted_auc_network2 <- targeted_removal_random_ties(network2, iterations = 1000)
targeted_auc_network3 <- targeted_removal_random_ties(network3, iterations = 1000)

num_iterations <- 1000 #declare for dataset combining
# Combine iteration data with mean values for the barplot
iteration_data <- data.frame(
  Network = factor(rep(c("Network 1", "Network 2", "Network 3"), each = 2 * num_iterations)),
  Strategy = factor(rep(c("Random", "Targeted (Ties)"), each = num_iterations, times = 3)),
  AUC = c(random_auc_network1, targeted_auc_network1,
          random_auc_network2, targeted_auc_network2,
          random_auc_network3, targeted_auc_network3)
)

# Summary data for bars (mean AUC values)
auc_comparison <- data.frame(
  Network = factor(rep(c("Network 1", "Network 2", "Network 3"), each = 2)),
  Strategy = factor(rep(c("Random", "Targeted (Ties)"), times = 3)),
  Mean_AUC = c(mean(random_auc_network1), mean(targeted_auc_network1),
               mean(random_auc_network2), mean(targeted_auc_network2),
               mean(random_auc_network3), mean(targeted_auc_network3)),
  SD = c(sd(random_auc_network1), sd(targeted_auc_network1),
         sd(random_auc_network2), sd(targeted_auc_network2),
         sd(random_auc_network3), sd(targeted_auc_network3)) # Standard deviation
)
auc_comparison



library(ggplot2)

# Filter data for Random and Targeted strategies
iteration_data_random <- subset(iteration_data, Strategy == "Random")
iteration_data_targeted <- subset(iteration_data, Strategy == "Targeted (Ties)")

auc_comparison_random <- subset(auc_comparison, Strategy == "Random")
auc_comparison_targeted <- subset(auc_comparison, Strategy == "Targeted (Ties)")

# Random Removal Plot
plot_random <- ggplot() +
  # Barplot for mean AUC (Random)
  geom_bar(data = auc_comparison_random, aes(x = Network, y = Mean_AUC, fill = Strategy),
           stat = "identity", width = 0.7, color = "black") +
  # Overlay iteration points (Random)
  geom_jitter(data = iteration_data_random, aes(x = Network, y = AUC),
              color = "darkblue", position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  # Add error bars
  geom_errorbar(data = auc_comparison_random, aes(x = Network, ymin = Mean_AUC - SD, ymax = Mean_AUC + SD),
                width = 0.25) +
  # Set y-axis limits
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = "AUC Comparison: Random Removal",
    x = "Network",
    y = "AUC"
  ) +
  theme_classic() +
  scale_fill_manual(values = "steelblue") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
plot_random


# Targeted Removal Plot
plot_targeted <- ggplot() +
  # Barplot for mean AUC (Targeted)
  geom_bar(data = auc_comparison_targeted, aes(x = Network, y = Mean_AUC, fill = Strategy),
           stat = "identity", width = 0.7, color = "black") +
  # Overlay iteration points (Targeted)
  geom_jitter(data = iteration_data_targeted, aes(x = Network, y = AUC),
              color = "darkred", position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
  # Add error bars
  geom_errorbar(data = auc_comparison_targeted, aes(x = Network, ymin = Mean_AUC - SD, ymax = Mean_AUC + SD),
                width = 0.25) +
  # Set y-axis limits
  scale_y_continuous(limits = c(0, 0.8)) +
  labs(
    title = "AUC Comparison: Targeted Removal with Tie-Breaking",
    x = "Network",
    y = "AUC"
  ) +
  theme_classic() +
  scale_fill_manual(values = "tomato") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

plot_targeted


library(car) # For Levene's test
library(multcomp) # For pairwise comparisons (Tukey's HSD)

# Function to perform tests and post-hoc analyses
perform_analysis <- function(data, strategy_label) {
  print(paste("### Analysis for", strategy_label, "###"))
  
  # Test for normality
  print("Shapiro-Wilk Test for Normality:")
  shapiro_results <- by(data$AUC, data$Network, shapiro.test)
  print(shapiro_results)
  
  # Check if normality holds
  normality_holds <- all(sapply(shapiro_results, function(x) x$p.value > 0.05))
  
  if (normality_holds) {
    # Levene's Test for Equal Variances
    print("Levene's Test for Equal Variances:")
    levene_result <- leveneTest(AUC ~ Network, data = data)
    print(levene_result)
    
    if (levene_result$`Pr(>F)`[1] > 0.05) {
      # Perform ANOVA
      print("Performing ANOVA:")
      anova_result <- aov(AUC ~ Network, data = data)
      print(summary(anova_result))
      
      # Post-hoc test: Tukey's HSD
      print("Tukey's HSD Post-Hoc Test:")
      tukey_result <- TukeyHSD(anova_result)
      print(tukey_result)
      
      # Extract adjusted p-values
      tukey_pvalues <- as.data.frame(tukey_result$Network)
      tukey_pvalues$Comparison <- rownames(tukey_result$Network)
      print(tukey_pvalues)
    } else {
      print("Variances are not equal; switching to Kruskal-Wallis test.")
      # Perform Kruskal-Wallis Test
      kruskal_result <- kruskal.test(AUC ~ Network, data = data)
      print(kruskal_result)
      
      # Post-hoc test: Pairwise Wilcoxon
      print("Pairwise Wilcoxon Post-Hoc Test with Bonferroni Correction:")
      pairwise_result <- pairwise.wilcox.test(data$AUC, data$Network, p.adjust.method = "fdr")
      print(pairwise_result)
    }
  } else {
    # Perform Kruskal-Wallis Test
    print("Performing Kruskal-Wallis Test:")
    kruskal_result <- kruskal.test(AUC ~ Network, data = data)
    print(kruskal_result)
    
    # Post-hoc test: Pairwise Wilcoxon
    print("Pairwise Wilcoxon Post-Hoc Test with Bonferroni Correction:")
    pairwise_result <- pairwise.wilcox.test(data$AUC, data$Network, p.adjust.method = "fdr")
    print(pairwise_result)
  }
}

# Perform analysis for Random and Targeted strategies
perform_analysis(iteration_data_random, "Random Removal") 
perform_analysis(iteration_data_targeted, "Targeted Removal with Tie-Breaking")




