# author:
#   Adam Gordon-Fennell (agg2248@uw.edu)
#   2022
#   Garret Stuber Lab, University of Washington

import sys
sys.path.append('./functions')
from py_fp import tidy_tdt_extract_and_tidy


# suggested format
#  - windows: r'C:\...\final_folder'
#  - mac: '/Users/.../final_folder/'

dir_raw = r'C:\Users\stuberadmin\Documents\GitHub\tidy_lab_tools_python\examples\tdt\raw'
dir_extracted = r'C:\Users\stuberadmin\Documents\GitHub\tidy_lab_tools_python\examples\tdt\extracted'


tidy_tdt_extract_and_tidy(dir_raw, dir_extracted)
