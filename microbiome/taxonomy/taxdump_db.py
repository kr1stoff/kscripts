#!/usr/bin/env python
import os
import sys
import sqlite3
import tarfile
import shutil
import pdb


def create_taxdump_db(db_dir, taxdump_tar):
    """
    Create taxdump sqlite3 db
    input: 
        [0] output sqlite db path; [1] directory which uncompressed taxdump.tar.gz 
    output:
        taxdump sqlite3 db
        >>>
        sqlite> SELECT * FROM taxdump LIMIT 5;
        tax_id  name                      rank          parent_tax_id  old_tax_id
        ------  ------------------------  ------------  -------------  ----------
        1       root                      no rank       1              NA        
        2       Bacteria                  superkingdom  131567         NA        
        6       Azorhizobium              genus         335928         NA        
        7       Azorhizobium caulinodans  species       6              395       
        9       Buchnera aphidicola       species       32199          28241 
        <<<
    """
    taxdump_tar_tmpdir = os.path.join(db_dir, "taxdump_tmp")
    #~ if taxdump_tar_tmpdir exists, remove it.
    if os.path.exists(taxdump_tar_tmpdir):
        shutil.rmtree(taxdump_tar_tmpdir)
    #~ make temporary taxdump directory, decompression 
    os.mkdir(taxdump_tar_tmpdir)
    with tarfile.open(taxdump_tar, "r:gz") as tar:
        tar.extractall(path=taxdump_tar_tmpdir)
    #~ create nodes names merged dictionary
    names_dict, nodes_dict, merged_dict = create_nodes_names_merged_dict(taxdump_tar_tmpdir=taxdump_tar_tmpdir)
    #~ sqlite3
    taxdump_db_file = os.path.join(db_dir, "taxdump.db")
    if os.path.exists(taxdump_db_file):
        os.remove(taxdump_db_file)
    con = sqlite3.connect(taxdump_db_file)
    cur = con.cursor()
    cur.execute('''CREATE TABLE taxdump
                (tax_id text, 
                name text, 
                rank text, 
                parent_tax_id text, 
                old_tax_id text)''')
    for tax_id in nodes_dict.keys():
        parent_tax_id, rank = nodes_dict[tax_id]
        # if-else one-line assign variate
        name = names_dict[tax_id]           #if tax_id in names_dict else "NA"
        old_tax_id = merged_dict[tax_id]    if tax_id in merged_dict else "NA"
        cur.execute("INSERT INTO taxdump VALUES (?, ?, ?, ?, ?)", (tax_id, name, rank, parent_tax_id, old_tax_id))
    con.commit()
    con.close()

