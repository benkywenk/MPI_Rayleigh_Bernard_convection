#Merge distributed analysis sets from a FileHandler.
import pathlib
from dedalus.tools import post

import logging

logger = logging.getLogger(__name__)

logger.info('Merging data files')

#Merging snapshot files
post.merge_process_files(save_direc + "snapshots", cleanup=True)
set_paths = list(pathlib.Path(save_direc + "snapshots/").glob("snapshots_s*.h5"))
post.merge_sets(save_direc + "snapshots/snapshots.h5", set_paths, cleanup=True)

#Merging analysis files
post.merge_process_files(save_direc + "analysis", cleanup=True)
set_paths = list(pathlib.Path(save_direc + "analysis/").glob("analysis_s*.h5"))
post.merge_sets(save_direc + "analysis/analysis.h5", set_paths, cleanup=True)

#Merging run_parameter files
post.merge_process_files(save_direc + "run_parameters", cleanup=True)
set_paths = list(pathlib.Path(save_direc + "run_parameters/").glob("run_parameters_s*.h5"))
post.merge_sets(save_direc + "run_parameters/run_parameters.h5", set_paths, cleanup=True)
