#!/usr/bin/env python
import os
import re
from collections import defaultdict
from glob import glob

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from ete3 import BarChartFace, TextFace, Tree, TreeStyle
from ete3.treeview.faces import add_face_to_node

os.environ["QT_QPA_PLATFORM"] = "offscreen"


def add_support_values(tree):
    for node in tree.traverse():
        if not node.is_leaf() and "/" in node.name:
            shalrt, bootstrap = [float(x) for x in node.name.split("/")]
            node.add_features(
                SHaLRT=shalrt,
                bootstrap=bootstrap,
                barchart_values=defaultdict(lambda: 0),
            )
            node.name = ""


outgroup = "K.CD.87.P3844.MH705156"
bootstrap_cutoff = 95

whole_aln_tree = Tree("trees/HIV1_FLT_2018_genome_DNA.fa.treefile", format=1)
whole_aln_tree.set_outgroup(outgroup)
add_support_values(whole_aln_tree)

mask_regex = r"mask(\d+)"
pct_masks = [0]
node_support_values = []
tree_file_list = [t for t in glob("trees/HIV1_FLT_2018_genome_DNA_mask*.fa.treefile")]
for tree_file in tree_file_list:
    pct_mask_match = re.search(r"mask(\d+)", tree_file)
    if pct_mask_match is not None:
        print("Adding {} bootstrap values to whole alignment tree...".format(tree_file))
        pct_mask = int(pct_mask_match.groups()[0])
        pct_masks.append(pct_mask)
        masked_tree = Tree(tree_file, format=1)
        add_support_values(masked_tree)
        masked_tree.set_outgroup(outgroup)
        comparison = whole_aln_tree.compare(masked_tree)
        for common_edge in comparison["common_edges"]:
            masked_tree_node = masked_tree.get_common_ancestor(common_edge)
            if hasattr(masked_tree_node, "bootstrap"):
                whole_aln_tree_node = whole_aln_tree.get_common_ancestor(common_edge)
                whole_aln_tree_node.barchart_values[
                    pct_mask
                ] = masked_tree_node.bootstrap
                whole_aln_tree_node.barchart_values[0] = whole_aln_tree_node.bootstrap
                node_support_values.append(
                    {
                        "Whole Alignment Bootstrap": whole_aln_tree_node.bootstrap,
                        "Masked Bootstrap": masked_tree_node.bootstrap,
                        "Percent Mask": pct_mask,
                    }
                )
        print("Done.")


# plot trees


def botstrap_stars(node):
    if hasattr(node, "barchart_values") and len(node.barchart_values) > 1:
        if node.barchart_values[0] < bootstrap_cutoff and any(
            [node.barchart_values[x] > bootstrap_cutoff for x in pct_masks[1:]]
        ):
            add_face_to_node(
                face=TextFace("***"),
                node=node,
                column=0,
            )


def botstrap_bars(node):
    botstrap_stars(node)
    if hasattr(node, "barchart_values") and len(node.barchart_values) > 1:
        add_face_to_node(
            face=BarChartFace(
                [node.barchart_values[pct] for pct in sorted(pct_masks)],
                colors=[
                    "#34073D",
                    "#471240",
                    "#591D43",
                    "#6C2846",
                    "#7F3349",
                    "#923E4D",
                    "#A44850",
                    "#B75353",
                    "#CA5E56",
                    "#DC6959",
                    "#EF745C",
                ],
                max_value=100,
            ),
            node=node,
            column=0,
        )


stars_style = TreeStyle()
stars_style.layout_fn = botstrap_stars
whole_aln_tree.render(file_name="plots/tree_botstrap_stars.pdf", tree_style=stars_style)
bars_style = TreeStyle()
bars_style.layout_fn = botstrap_bars
whole_aln_tree.render(file_name="plots/tree_botstrap_bars.pdf", tree_style=bars_style)


node_support_df = pd.DataFrame(node_support_values)
# "top-right" of scatter plot
high_support_df = node_support_df[(node_support_df["Whole Alignment Bootstrap"] >=bootstrap_cutoff) & (node_support_df["Masked Bootstrap"]>=bootstrap_cutoff)]
whole_aln_support_df = node_support_df[(node_support_df["Whole Alignment Bootstrap"] >=bootstrap_cutoff) & (node_support_df["Masked Bootstrap"]<bootstrap_cutoff)]
masked_aln_support_df = node_support_df[(node_support_df["Whole Alignment Bootstrap"] <bootstrap_cutoff) & (node_support_df["Masked Bootstrap"]>=bootstrap_cutoff)]
low_support_df = node_support_df[(node_support_df["Whole Alignment Bootstrap"] <bootstrap_cutoff) & (node_support_df["Masked Bootstrap"]<bootstrap_cutoff)]

# scatter plot
color_map = mpl.cm.get_cmap('cividis')
#.colors.to_hex(
#color_map = mpl.cm.get_cmap('cubehelix')
scatter_pcts = node_support_df["Percent Mask"].unique()
scatter_pcts.sort()
scatter_colors = color_map((scatter_pcts)/100)
fig, ax = plt.subplots(dpi=300)
ax.set_facecolor("gainsboro")
xlabel = "Whole Alignment Bootstrap"
ylabel = "Masked Bootstrap"
for pct, color in zip(scatter_pcts[::-1], scatter_colors[::-1]):
    pct_df = node_support_df[node_support_df["Percent Mask"]==pct]
    ax.scatter(pct_df[xlabel], pct_df[ylabel], color=color, label="{}% masked".format(pct), alpha=0.4, edgecolors='none')
line_styles = {"color":'black', "linestyle":"--"}
ax.axvline(bootstrap_cutoff, **line_styles)
ax.axhline(bootstrap_cutoff, **line_styles)

top_right_pct_nodes = high_support_df.shape[0]/node_support_df.shape[0]*100
bottom_right_pct_nodes = whole_aln_support_df.shape[0]/node_support_df.shape[0]*100
top_left_pct_nodes = masked_aln_support_df.shape[0]/node_support_df.shape[0]*100
bottom_left_pct_nodes = low_support_df.shape[0]/node_support_df.shape[0]*100
ax.text(bootstrap_cutoff+2, bootstrap_cutoff+2, "{:.1f}%".format(top_right_pct_nodes))
ax.text(bootstrap_cutoff+2, 2, "{:.1f}%".format(bottom_right_pct_nodes))
ax.text(2, bootstrap_cutoff+2, "{:.1f}%".format(top_left_pct_nodes))
ax.text(2,2,"{:.1f}%".format(bottom_left_pct_nodes))

ax.set_xlim([0, 102])
ax.set_xlabel(xlabel)
ax.set_ylim([0, 102])
ax.set_ylabel(ylabel)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig('plots/bootstrap_scatter.png')
