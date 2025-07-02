# author:
#   Adam Gordon-Fennell (agg2248@uw.edu)
#   2022
#   Garret Stuber Lab, University of Washington

import sys
sys.path.append('./functions')
from py_fp import tidy_doric_extract_and_tidy


dir_raw = r'C:\Users\stuberadmin\Desktop\temp\doric\batch\raw'
dir_extracted = r'C:\Users\stuberadmin\Desktop\temp\doric\batch\extracted'

tidy_doric_extract_and_tidy(dir_raw, dir_extracted)

