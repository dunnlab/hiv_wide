#!/usr/bin/env python
import re
from glob import glob
from collections import defaultdict

from ete3 import BarChartFace, TextFace, Tree, TreeStyle
from ete3.treeview.faces import add_face_to_node


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
# common_bootstraped_nodes = set()
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
                whole_aln_tree_node.barchart_values[0] = masked_tree_node.bootstrap
                # common_bootstraped_nodes.add(whole_aln_tree_node)
        print("Done.")


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
                #min_value=80,
                max_value=100,
            ),
            node=node,
            column=0,
        )



stars_style = TreeStyle()
stars_style.layout_fn = botstrap_stars
whole_aln_tree.render(file_name="trees/botstrap_stars.pdf", tree_style=stars_style)
bars_style = TreeStyle()
bars_style.layout_fn = botstrap_bars
whole_aln_tree.render(file_name="trees/botstrap_bars.pdf", tree_style=bars_style)
