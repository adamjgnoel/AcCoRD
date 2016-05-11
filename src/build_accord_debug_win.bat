@echo off
if not exist "..\bin\" mkdir "..\bin"
gcc accord.c pcg_basic.c err_fcts.c w_of_z.c im_w_of_x.c erfcx.c rand_accord.c region.c subvolume.c meso.c micro_molecule.c chem_rxn.c actor.c base.c observations.c timer_accord.c mol_release.c actor_data.c file_io.c cJSON.c -std=c99 -g -o "..\bin\accord_win_debug.exe"