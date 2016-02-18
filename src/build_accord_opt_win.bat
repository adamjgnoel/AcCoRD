@echo off
if not exist "..\bin\" mkdir "..\bin"
gcc accord.c mtwist.c randistrs.c subvolume.c meso.c micro_molecule.c region.c chem_rxn.c actor.c base.c observations.c timer_accord.c mol_release.c actor_data.c file_io.c cJSON.c -std=c99 -O3 -o "..\bin\accord_win.exe"