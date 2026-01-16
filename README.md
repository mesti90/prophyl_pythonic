# prophyl_pythonic
Basically the same pipeline as kintses_lab/prophyl, but it is built in python instead of nextflow

**pythonic_prophyl.py**
*Necessary args:*
- st_file: the list of STs to be processed
- typing_table: a tab-sep table containing all metadata. Two important columns: strain and Date. For the figures, we also use spatyper, Capsule type and country columns
- all_assemblies: a tab-sep table containing the paths of the assemblies. Columns used: strain and assembly_path
