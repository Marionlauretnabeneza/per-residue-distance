
# Plotting code for per-residue atomic distances



Obtain rotation-translation matrices from the PDBe FTP server:

```https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition/<uniprot>[0]/<uniprot>/<uniprot>.json```

E.g. 

```wget https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition/A/A0QTT2/A0QTT2.json```


Obtain updated mmCIF files (the normal, non-SIFTS mmCIF file type will not work) using:


```https://www.ebi.ac.uk/pdbe/entry-files/download/<pdbid>_updated.cif```

E.g.

```wget https://www.ebi.ac.uk/pdbe/entry-files/download/7cyr_updated.cif```

## Example 

``` python
python3 per_residue_distance.py --mmcif1 example_data/7cyr_updated.cif.gz --mmcif2 example_data/7cy2_updated.cif.gz --rt_matrices example_data/A0QTT2.json --pdb_id1 7cyr --pdb_id2 7cy2 --chain1 A --chain2 A
```

Chain IDs should be parsed as `label_asym_id` from the updated mmcif `atom_site` loop. 