#!/usr/bin/env python
import os
import re
import sys
from collections import defaultdict
from glob import glob

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from ete3 import BarChartFace, NodeStyle, RectFace, TextFace, Tree, TreeStyle
from ete3.treeview.faces import add_face_to_node

os.environ["QT_QPA_PLATFORM"] = "offscreen"
outgroup = "K.CD.87.P3844.MH705156"
bootstrap_cutoff = 95


def add_support_and_subtypes(tree):
    for node in tree.traverse():
        if not node.is_leaf() and "/" in node.name:
            shalrt, bootstrap = [float(x) for x in node.name.split("/")]
            node.add_features(
                SHaLRT=shalrt,
                bootstrap=bootstrap,
                barchart_values=defaultdict(lambda: 0),
                suport_symbol="",
            )
            node.name = ""
        elif node.is_leaf():
            subtype = node.name.split('.',1)[0]
            node.add_features(
                subtype = subtype
            )


def get_support_symbol(ref_aln_bstrap, comp_aln_bstrap):
    if ref_aln_bstrap >= bootstrap_cutoff and comp_aln_bstrap >= bootstrap_cutoff:
        return "↗"
    elif ref_aln_bstrap >= bootstrap_cutoff and comp_aln_bstrap < bootstrap_cutoff:
        return "↘"
    elif ref_aln_bstrap < bootstrap_cutoff and comp_aln_bstrap >= bootstrap_cutoff:
        return "↖"
    elif ref_aln_bstrap < bootstrap_cutoff and comp_aln_bstrap < bootstrap_cutoff:
        return "↙"


ref_tree_file = "trees/HIV1_FLT_2018_genome_DNA.fa.treefile"
plot_prefix = "whole_v_masked"
ref_bs_label = "Whole Alignment Bootstrap"
tree_bs_label = "Masked Alignment Bootstrap"
tree_file_list = ["trees/HIV1_FLT_2018_genome_DNA_mask100.fa.treefile"]
tree_orientation = 0
pct_masks = [0]

if "--all" in sys.argv:
    tree_file_list = [
        t for t in glob("trees/HIV1_FLT_2018_genome_DNA_mask*.fa.treefile")
    ]
    plot_prefix = "whole_v_allmasks"
if "--mask-as-ref" in sys.argv:
    tree_file_list = [ref_tree_file]
    ref_tree_file = "trees/HIV1_FLT_2018_genome_DNA_mask100.fa.treefile"
    plot_prefix = "masked_v_whole"
    ref_bs_label = "Masked Alignment Bootstrap"
    tree_bs_label = "Whole Alignment Bootstrap"
    pct_masks = [100]
    tree_orientation = 1


ref_tree = Tree(ref_tree_file, format=1)
ref_tree.set_outgroup(outgroup)
add_support_and_subtypes(ref_tree)

mask_regex = r"mask(\d+)"

shared_edge_support_values = []
ref_only_edge_support_values = []
for tree_file in tree_file_list:
    pct_mask_match = re.search(r"mask(\d+)", tree_file)
    if pct_mask_match is not None:
        pct_mask = int(pct_mask_match.groups()[0])
    else:
        pct_mask = 0
    print("Adding {} bootstrap values to tree from {}...".format(tree_file, ref_tree_file))
    pct_masks.append(pct_mask)
    tree = Tree(tree_file, format=1)
    add_support_and_subtypes(tree)
    tree.set_outgroup(outgroup)
    comparison = ref_tree.compare(tree)
    print("{}/{} common/total edges, normRF {:0.2f} for {} vs {}".format(len(comparison["common_edges"]), len(comparison["source_edges"]), comparison["norm_rf"], tree_file, ref_tree_file ))
    for common_edge in comparison["common_edges"]:
        tree_node = tree.get_common_ancestor(common_edge)
        if hasattr(tree_node, "bootstrap"):
            ref_tree_node = ref_tree.get_common_ancestor(common_edge)
            ref_tree_node.barchart_values[pct_mask] = tree_node.bootstrap
            ref_tree_node.barchart_values[0] = ref_tree_node.bootstrap
            ref_tree_node.suport_symbol = get_support_symbol(
                ref_tree_node.bootstrap, tree_node.bootstrap
            )
            shared_edge_support_values.append(
                {
                    ref_bs_label: ref_tree_node.bootstrap,
                    tree_bs_label: tree_node.bootstrap,
                    "Percent Mask": pct_mask,
                }
            )
    for ref_only_edge in comparison["ref_edges"] - comparison["common_edges"]:
        node = ref_tree.get_common_ancestor(ref_only_edge)
        if hasattr(node, "bootstrap"):
            ref_only_edge_support_values.append(node.bootstrap)

    print("Done.")


# Make some colors

subtypes = ["A", "B", "C", "D", "F1", "F2", "G", "H", "J", "K"]
tab10_cmap = mpl.cm.get_cmap("tab10")
subtype_color_dict = dict(zip(subtypes, [mpl.colors.to_hex(x) for x in tab10_cmap(np.linspace(0,1,10))] ))

color_map = mpl.cm.get_cmap("cividis")
plot_colors = color_map([x / 100 for x in reversed(sorted(pct_masks))])
tree_colors = [mpl.colors.to_hex(x) for x in plot_colors]
scatter_colors = plot_colors

# plot trees
# color monophyletic subtypes
print("Plotting trees...")
def color_subtypes(node):
    node_style = NodeStyle()
    node_style["hz_line_width"] = 2
    node_style["vt_line_width"] = 2
    if hasattr(node, "bootstrap"):
        if node.bootstrap >= bootstrap_cutoff:
            node_style["fgcolor"] = "black"

        else:
            node_style["fgcolor"] = "grey"
    else:
        node_style["fgcolor"] = "black"
    
    if hasattr(node, "subtype") and node.subtype in subtype_color_dict:
        node_style["hz_line_color"] = subtype_color_dict[node.subtype]
    else:
        # check for subtype monophyly and color if true
        subtypes = set()
        for leaf in node.get_leaves():
            subtypes.add(leaf.subtype)
            if len(subtypes) > 1:
                break
        subtype = subtypes.pop()
        if len(subtypes) == 0 and subtype in subtype_color_dict:
            node_style["vt_line_color"] = subtype_color_dict[subtype]
            node_style["hz_line_color"] = subtype_color_dict[subtype]
    return node_style

# add an arrow symbol pointing to the quadrant 
# this bootstrap comparison belongs in in scatter below
def botstrap_symbols(node):
    node_style = color_subtypes(node)
    if hasattr(node, "suport_symbol"):
        add_face_to_node(
            face=TextFace(node.suport_symbol),
            node=node,
            column=0,
        )
    node.set_style(node_style)

# call larger attention to lower right quadrant from scatter
# 
def botstrap_lower_right(node):
    node_style = color_subtypes(node)
    if hasattr(node, "suport_symbol"):
        if node.suport_symbol == "↘":
            node_style["size"] = 100
            node_style["fgcolor"] = "#66000000"

    node.set_style(node_style)


for output_string, mode, layout_fn in zip(["tree_symbols", "tree_symbols", "tree_arrows",], ["c", "r", "r"], [botstrap_lower_right, botstrap_lower_right, botstrap_symbols]):
    stars_style = TreeStyle()
    stars_style.layout_fn = layout_fn
    stars_style.orientation = tree_orientation
    stars_style.mode = mode
    for subtype in subtype_color_dict:
        stars_style.legend.add_face(RectFace(10,10, subtype_color_dict[subtype], subtype_color_dict[subtype]),column=0)
        stars_style.legend.add_face(TextFace(subtype, fgcolor=subtype_color_dict[subtype], bold=True),column=1)
    stars_style.legend_position = 2
    ref_tree.render(
        file_name="plots/{}_{}_{}.pdf".format(plot_prefix, output_string, mode), tree_style=stars_style
    )
print("Done.")

# dataframe for bootstrap values plotting
shared_edge_df = pd.DataFrame(shared_edge_support_values)
# "top-right" of scatter plot
high_support_df = shared_edge_df[
    (shared_edge_df[ref_bs_label] >= bootstrap_cutoff)
    & (shared_edge_df[tree_bs_label] >= bootstrap_cutoff)
]
# "bottom-right" of scatter plot
ref_aln_support_df = shared_edge_df[
    (shared_edge_df[ref_bs_label] >= bootstrap_cutoff)
    & (shared_edge_df[tree_bs_label] < bootstrap_cutoff)
]
# "top-left" of scatter plot
other_aln_support_df = shared_edge_df[
    (shared_edge_df[ref_bs_label] < bootstrap_cutoff)
    & (shared_edge_df[tree_bs_label] >= bootstrap_cutoff)
]
# "bottom-left" of scatter plot
low_support_df = shared_edge_df[
    (shared_edge_df[ref_bs_label] < bootstrap_cutoff)
    & (shared_edge_df[tree_bs_label] < bootstrap_cutoff)
]

# scatter plot
print("Plotting scatter...")
fig, ax = plt.subplots(dpi=300)
ax.set_facecolor("gainsboro")
xlabel = ref_bs_label
ylabel = tree_bs_label
scatter_pcts = shared_edge_df["Percent Mask"].unique()
scatter_pcts.sort()
for pct, color in zip(reversed(scatter_pcts), reversed(scatter_colors)):
    pct_df = shared_edge_df[shared_edge_df["Percent Mask"] == pct]
    ax.scatter(
        pct_df[xlabel],
        pct_df[ylabel],
        color=color,
        label="{}% masked".format(pct),
        alpha=0.4,
        edgecolors="none",
    )
line_styles = {"color": "black", "linestyle": "--"}
ax.axvline(bootstrap_cutoff, **line_styles)
ax.axhline(bootstrap_cutoff, **line_styles)

top_right_pct_nodes = high_support_df.shape[0] / shared_edge_df.shape[0] * 100
bottom_right_pct_nodes = ref_aln_support_df.shape[0] / shared_edge_df.shape[0] * 100
top_left_pct_nodes = other_aln_support_df.shape[0] / shared_edge_df.shape[0] * 100
bottom_left_pct_nodes = low_support_df.shape[0] / shared_edge_df.shape[0] * 100
ax.text(
    101, bootstrap_cutoff + 2, "{:.1f}%".format(top_right_pct_nodes)
)
ax.text(bootstrap_cutoff + 2, 2, "{:.1f}%".format(bottom_right_pct_nodes))
ax.text(2, bootstrap_cutoff + 2, "{:.1f}%".format(top_left_pct_nodes))
ax.text(2, 2, "{:.1f}%".format(bottom_left_pct_nodes))

ax.set_xlim([0, 102])
ax.set_xlabel(xlabel)
ax.set_ylim([0, 102])
ax.set_ylabel(ylabel)
if len(scatter_pcts) >1:
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
plt.tight_layout()
plt.savefig("plots/{}_bootstrap_scatter.png".format(plot_prefix))
plt.close()
print("Done.")

print("Plotting histogram...")
fig, ax = plt.subplots(dpi=300)
ax.hist(ref_only_edge_support_values, bins=np.arange(10,105,5))
ax.set_xlabel("{}".format(xlabel))
plt.tight_layout()
plt.savefig("plots/{}_bootstrap_hist.png".format(plot_prefix))
print("Done.")

"""
# quick whole/pol hist
%run ./scripts/bootstrap_support.py
whole_aln = ref_only_edge_support_values
%run ./scripts/bootstrap_support.py --mask-as-ref
masked_aln = ref_only_edge_support_values
plt.close()

fig, ax = plt.subplots(dpi=300)
ax.hist([masked_aln, whole_aln], bins=np.arange(10,105,5), label=["pol", "whole"])
ax.set_xlabel("Bootstrap")
plt.legend()
plt.tight_layout()
plt.savefig("plots/wholevpol_hist.png")
"""