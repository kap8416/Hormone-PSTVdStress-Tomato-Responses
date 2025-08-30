# Load required libraries
library(igraph)
library(readr)
library(dplyr)

# Set seed for reproducibility
set.seed(22)

# Define base colors
color_ON <- "#91ee91"   # Color for ON nodes
color_OFF <- "#A3195B"  # Color for OFF nodes

# Design parameters
default_border_on <- "#865bb3"  # Border for ON nodes
default_border_off <- "#6A6C70" # Border for OFF nodes

# Function to create graph with mode-based coloring
create_static_network <- function(df_edges,
                                  border_on = default_border_on,
                                  border_off = default_border_off) {
  
  # Create node dataframe
  nodes_df <- unique(data.frame(
    node = c(df_edges$node1, df_edges$node2),
    mode = c(df_edges$node1_mode, df_edges$node2_mode),
    stringsAsFactors = FALSE
  ))
  
  # Assign fill colors based on mode only
  nodes_df$color <- ifelse(nodes_df$mode == "1", color_ON, color_OFF)
  
  # Create graph
  g <- graph_from_data_frame(
    d = df_edges[, c("node1", "node2")],
    directed = TRUE,
    vertices = nodes_df
  )
  
  # Assign border colors based on mode only
  V(g)$frame.color <- ifelse(V(g)$mode == "1", border_on, border_off)
  
  # Set border width
  V(g)$frame.width <- 2
  
  # Edge color logic
  E(g)$color <- apply(get.edges(g, 1:ecount(g)), 1, function(edge) {
    source <- V(g)[edge[1]]
    target <- V(g)[edge[2]]
    
    if (source$mode == "1" && target$mode == "1") {
      if (target$name %in% df_edges$node1) "red" else "#6A6C70"
    } else {
      "gray"
    }
  })
  
  # Edge width logic
  E(g)$width <- ifelse(E(g)$color == "#6A6C70", 3, 
                       ifelse(E(g)$color == "red", 2, 1))
  
  return(g)
}

# Function to plot and save network
plot_network <- function(g, layout, file_prefix, show_labels = TRUE) {
  # Get edges and vertices
  edges_df <- get.data.frame(g, what = "edges")
  vertices_df <- get.data.frame(g, what = "vertices")
  
  # Reorder edges (red on top)
  edges_df$edge_order <- ifelse(edges_df$color == "red", 2, 1)
  edges_df <- edges_df[order(edges_df$edge_order), ]
  
  # Rebuild graph with reordered edges
  g_new <- graph_from_data_frame(edges_df, directed = TRUE, vertices = vertices_df)
  
  # Save PNG with transparent background
  png(paste0(file_prefix, ".png"), width = 2000, height = 1500, res = 300, bg = "transparent")
  plot(g_new, layout = layout,
       vertex.size = 15,
       vertex.label = if (show_labels) V(g_new)$name else NA,
       vertex.label.cex = 0.8,
       edge.arrow.size = 0.5,
       vertex.frame.width = V(g_new)$frame.width,
       background = NA)
  dev.off()
  
  # Save SVG with transparent background
  svg(paste0(file_prefix, ".svg"), width = 20, height = 15, bg = "transparent")
  plot(g_new, layout = layout,
       vertex.size = 15,
       vertex.label = if (show_labels) V(g_new)$name else NA,
       vertex.label.cex = 0.8,
       edge.arrow.size = 0.5 * 5,
       edge.width = E(g_new)$width * 2.5,
       vertex.frame.width = V(g_new)$frame.width * 6,
       background = NA)
  dev.off()
}

# Process all TSV files in folder
folder_path <- "./"  # Change to your directory
tsv_files <- list.files(path = folder_path, pattern = "*.tsv", full.names = TRUE)

# List to store graphs
graph_list <- list()

for (file in tsv_files) {
  # Read TSV without headers (4 columns)
  df_edges <- read_tsv(file, col_names = FALSE) %>%
    setNames(c("node1", "node2", "node1_mode", "node2_mode")) %>%
    mutate(across(ends_with("_mode"), as.character))
  
  # Create network
  g <- create_static_network(df_edges,
                             border_on = default_border_on,
                             border_off = default_border_off)
  
  # Add graph to list
  graph_list[[file]] <- g
}

# Calculate layout using first graph
base_graph <- graph_list[[1]]
layout_coordinates <- layout_with_fr(base_graph)

# Generate plots for all graphs
for (file_path in names(graph_list)) {
  g <- graph_list[[file_path]]
  base_name <- tools::file_path_sans_ext(basename(file_path))
  plot_network(g, layout_coordinates, base_name, show_labels = FALSE) #choose between labeled (TRUE) or unlabeled (FALSE) output networks
}