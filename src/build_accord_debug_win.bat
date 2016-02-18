@echo off
if not exist "..\bin\" mkdir "..\bin"
gcc accord.c mtwist.c randistrs.c region.c subvolume.c meso.c micro_molecule.c chem_rxn.c actor.c base.c observations.c timer_accord.c mol_release.c actor_data.c file_io.c cJSON.c -std=c99 -g -o "..\bin\accord_win_debug.exe"