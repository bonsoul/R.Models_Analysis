  # if not already installed
library(networkD3)
library(dplyr)

# Built-in dataset: Titanic (4-way contingency table)
df <- as.data.frame(Titanic)

# Build links: Class -> Sex -> Survived, weighted by Freq
links <- bind_rows(
  df %>% group_by(source = Class, target = Sex) %>% summarise(value = sum(Freq), .groups = "drop"),
  df %>% group_by(source = Sex, target = Survived) %>% summarise(value = sum(Freq), .groups = "drop")
)

links$source <- as.character(links$source)
links$target <- as.character(links$target)

# Node list
nodes <- data.frame(name = unique(c(links$source, links$target)))

# Convert names to zero-based indices for networkD3
links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", NodeID = "name",
  fontSize = 12, nodeWidth = 30
)
