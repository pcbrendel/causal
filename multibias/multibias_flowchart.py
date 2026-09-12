from graphviz import Digraph

# Create a directed graph
dot = Digraph(comment="Multibias Flow", format="png")

# Graph attributes for clean academic look
dot.attr(rankdir="TB", size="8", nodesep="0.5", ranksep="0.75")
dot.attr(
    "node",
    shape="rectangle",
    style="rounded,filled",
    fontname="Helvetica",
    fontsize="12",
    fillcolor="white",
)

# Nodes
dot.node(
  "data_obs",
  "data_observed",
  shape = "box",
  fillcolor="white"
)
dot.node(
  "bias",
  "bias_params\nor\ndata_validation",
  shape = "box",
  fillcolor="white"
)
dot.node(
  "mb_adjust",
  "1. multibias_adjust()",
  shape = "ellipse",
  fillcolor="lightgrey"
)
dot.node(
  "results",
  "Adjusted Estimate 1\nAdjusted Estimate 2\n...\nAdjusted Estimate n", fillcolor="white"
)
dot.node(
  "mb_plot",
  "2. multibias_plot()",
  shape = "ellipse",
  fillcolor="lightgrey"
)
dot.node(
  "note",
  "Run across different bias assumptions\n(bias parameters or validation data)",
  shape="note",
  fontsize="11",
  fillcolor="darkgrey",
)
dot.node(
  "forest",
  "Forest plot of adjusted estimates",
  shape="box",
  fillcolor="white",
)

# Edges
dot.edge("data_obs", "mb_adjust")
dot.edge("bias", "mb_adjust")
dot.edge("note", "results")
dot.edge("data_obs", "mb_plot")
dot.edge("results", "mb_plot")
dot.edge("mb_adjust", "note", style="dashed")
dot.edge("mb_plot", "forest")

# Save and render
output_path = "//Users/pbrendel/projects/causal/multibias_flow"
dot.render(output_path, cleanup=True)

output_path
