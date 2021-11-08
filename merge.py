#Merge distributed analysis sets from a FileHandler.
import pathlib
from dedalus.tools import post
import run_param_file as rpf

import logging

logger = logging.getLogger(__name__)

logger.info('Merging data files')
save_direc = "sim_data/"
run_name = rpf.run_name

#Merge function
def merge(dimensionality):
    #Merging snapshot files
    post.merge_process_files(save_direc + dimensionality + "snapshots", cleanup=True)
    set_paths = list(pathlib.Path(save_direc + dimensionality + "snapshots/").glob("snapshots_s*.h5"))
    post.merge_sets(save_direc + dimensionality + "snapshots/snapshots_" + run_name + ".h5", set_paths, cleanup=True)

    #Merging analysis files
    post.merge_process_files(save_direc + dimensionality + "analysis", cleanup=True)
    set_paths = list(pathlib.Path(save_direc + dimensionality + "analysis/").glob("analysis_s*.h5"))
    post.merge_sets(save_direc + dimensionality + "analysis/analysis_" + run_name + ".h5", set_paths, cleanup=True)

    #Merging run_parameter files
    post.merge_process_files(save_direc + dimensionality + "run_parameters", cleanup=True)
    set_paths = list(pathlib.Path(save_direc + dimensionality + "run_parameters/").glob("run_parameters_s*.h5"))
    post.merge_sets(save_direc + dimensionality + "run_parameters/run_parameters_" + run_name + ".h5", set_paths, cleanup=True)

#specifying simulation type
choice = input("Enter 3 for full 3D simulation 2 for 2.5D: ")

if choice == "3":
    merge("3D/")
elif choice == "2":
    merge("2.5D/")
else:
    print("Incorrect option, exiting program.")
