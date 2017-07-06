# Genome indexing using snakemake

Copy the ``Snakefile`` and ``config.json`` to a new directory. Then edit the ``config.json`` to reflect the desired genome build.
Execute Snakemake in the directory holding the newly created Snakefile like so:

```
  snakemake report -j 8
```

After the succesfull run you can inspect the content of
``add_{buildkey}_to_galaxy.sh``. Then as a Galaxy administrator you can execute
the script using bash to add all the indexes to your Galaxy installation.
