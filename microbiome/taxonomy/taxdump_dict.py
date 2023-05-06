#!/usr/bin/env python
import os
from collections import defaultdict


def dump2dict(taxdmp_dir:str, name_node=False):
    """
    Need uncompress taxdump.tar.gz to taxdmp_dir, change nodes names merged dump file to dictionaries.

    >>>
    taxdump_dict: {"nodes_dict": {..}, "names_dict": {..}, "merged_dict": {..}, ["name_node_dict": {..}]}

    nodes_dict: {
      "13": {
          "parent": "203488",
          "rank": "genus"}, 
    ...}
    names_dict: {"2": "Bacteria", ...}
    name_node_dict: {"Bacteria": ["2", ], ...} # some taxonomy hava same scientific name
    merged_dict: {"12": "74109", ...}
    <<<
    """
    nodes_dmp, names_dmp, merged_dmp = list(map(lambda x: os.path.join(taxdmp_dir, x), 
                                        ["nodes.dmp", "names.dmp", "merged.dmp"]))
    nodes_dict, names_dict, merged_dict, name_node_dict = dict(), dict(), dict(), defaultdict(list)
    #~ nodes
    with open(nodes_dmp, "rt") as f:
        for line in f:
            llist = line.replace("\t|\n", "").split("\t|\t")
            tax_id, parent_tax_id, rank = llist[:3]
            nodes_dict[tax_id] = {}
            nodes_dict[tax_id]["parent"] = parent_tax_id
            nodes_dict[tax_id]["rank"] = rank
    #~ names
    with open(names_dmp, "rt") as f:
        for line in f:
            #~ just scientific name
            if "scientific name" not in line:
                continue
            llist = line.replace("\t|\n", "").split("\t|\t")
            node, name = llist[:2]
            names_dict[node] = name
            name_node_dict[name].append(node)
    #~ merged
    with open(merged_dmp, "rt") as f:
        for line in f:
            llist = line.replace("\t|\n", "").split("\t|\t")
            old_tax_id, new_tax_id = llist
            merged_dict[old_tax_id] = new_tax_id
    #~ return output --> taxdump_dict: {"nodes_dict": {..}, "names_dict": {..}, "merged_dict": {..}, ["name_node_dict": {..}]}
    
    taxdump_dict = {
        "nodes_dict": nodes_dict,
        "names_dict": names_dict,
        "merged_dict": merged_dict
    }
    #` name - node dictionary, if need
    if name_node == True:
        taxdump_dict["name_node_dict"] = name_node_dict
    return taxdump_dict
