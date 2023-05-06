#!/usr/bin/env python


def exhibit_lineage(tax_id:str, nodes_dict:dict, merged_dict:dict, names_dict:dict):
    """
    Exhibit full lineage from kindom to current tax rank, such as kindom-species, kindom-family..."

    >>>
    output (USEARCH format):
        "d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus"
        or
        "no_record"
    <<<
    """
    lineage_list = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    lineage_enum_dict = {items[1]: items[0] for items in list(enumerate(lineage_list))}
    lineage_abbr_list = ["k","p","c","o","f","g","s"]
    #~ lineage_dict
    lineage_dict = dict()
    #~ check new_tax_id and old_tax_id, neither return "no_record"
    if tax_id in nodes_dict:
        tmp_tax_id = tax_id
    elif tax_id in merged_dict:
        tmp_tax_id = merged_dict[tax_id]
    else:
        return "no_record"
    #~ repetitive work 
    def add_for_lineage_dict(tax_id):
        tmp_rank = nodes_dict[tax_id]["rank"]
        tmp_name = names_dict[tax_id].replace(" ", "_")
        lineage_dict[tmp_rank] = tmp_name
    #~ main lineage_dict
    """
    {...
    'phylum': 'Actinobacteria',
    'clade': 'Terrabacteria group',
    'superkingdom': 'Bacteria',
    'no rank': 'cellular organisms'}
    """
    while True:
        add_for_lineage_dict(tmp_tax_id)
        tmp_tax_id = nodes_dict[tmp_tax_id]["parent"]
        if tmp_tax_id == "1":
            break
    #~ complete lineage_dict, set start rank (species, genus or higher..)
    tmp_rank_list = [rank for rank in lineage_list if rank in lineage_dict]
    last_rank_index = lineage_enum_dict[tmp_rank_list[-1]]
    final_rank_list = lineage_list[:(last_rank_index+1)]
    name_list = list()
    #` complete
    for rank in final_rank_list:
        if rank in lineage_dict:
            name_list.append(lineage_dict[rank])
        else:
            name_list.append("Undefined")
    final_lineage_abbr_list = lineage_abbr_list[:(last_rank_index+1)]
    #~ output full lineage, USEARCH format. Use commas to separate ranks, colons to name.
    # d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus
    out_full_lineage_list = list()
    for items in zip(final_lineage_abbr_list, name_list):
        out_full_lineage_list.append(":".join(items))
    out_full_lineage = ",".join(out_full_lineage_list)
    return out_full_lineage

# [disable] receive taxonomy id then find phynum
def find_phynum(tax_id, node_dict_with_rank, name_dict):
    """
    Disable! Fix the code if necessary.
    """
    flag = True
    #~ recursion  iteration
    tmp_tax_id = tax_id
    while flag:
        try:
            parent_tax_id = node_dict_with_rank[tmp_tax_id]['parent']
            rank = node_dict_with_rank[tmp_tax_id]['rank']
            if rank == 'phylum':
                flag = False
                phylum_name = name_dict[tmp_tax_id]
            elif parent_tax_id == '1':
                flag = False
            else:
                tmp_tax_id = parent_tax_id
        except NameError:
            flag = False
        except KeyError:
            flag = False
    if 'phylum_name' not in locals().keys():
        phylum_name = 'None'
    return phylum_name