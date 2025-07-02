## The structure of the code is as follows:
`functions/` contains both the R code and Python code. The Python code is separated into 2 files: `tdt_analysis.py` (where the bulk of the code is) and `general.py` (which contains a few helper functions).

The `tdt_analysis.py` file is an almost direct translation of the first 1260 lines of the `r_fp.R` file.

The code is outputted to `examples/tdt/processed` and reads from the `examples/tdt/extracted` folder and `examples/tdt/log_data_fp_tdt.csv` file

## File locations

* Example data: `\tidy_lab_tools_python\examples\tdt\extracted`
* Log_fp master list (user-defined list of sessions to analyze): `\tidy_lab_tools_python\examples\tdt`
* All-in-one GUI EXE (PREFERABLY USE THIS): `.\tidy_lab_tools_python\dist\gui\gui.exe`
* Main python code to run in IDE: `.\tidy_lab_tools_python\functions\tdt_analysis.py`
* Code to convert tdt raw data to tidy (data format required for GUI and tdt_analysis.py): `.\tidy_lab_tools_python\functions\tdt_to_tidy.py`

## Primary Steps

1. Convert any raw tdt or doric data to tidy format (tdt_to_tidy.py)
2. Place the data (named *_streams_data.csv or feather) and session meta info (*_streams_info.csv or feather) into the `.\tidy_lab_tools_python\examples\tdt\extracted` folder
3. Run Gui.exe or tdt_analysis.py; in the gui, navigate and select the `log_data_fp_tdt.csv` file
4. Visual inspection of traces
5. Use the output in `.\tidy_lab_tools_python\examples\tdt\processed` for downstream analyses

## Running the code

Easy approach:
Go to the `.\dist\gui\` folder and double click the gui.exe file

Python file approach:
Add data to the `examples/tdt/extracted` folder (if this folder does not exist, then create it), configure the log file `examples/tdt/log_data_fp_tdt_single_session.csv` (for example), and run `functions/tdt_analysis.py`

![image](https://github.com/user-attachments/assets/0a59ae95-ca62-4491-80aa-91df210ac3c4 =250x)

![image](https://github.com/user-attachments/assets/5760af11-8a65-4723-8b8b-8b7dc874337f)

## Data format:

data file (one for each session): 

* include columns with specific names: blockname, start_date, channel, and raw_au
* blockname, channel, and raw_au (raw fluorescence values) are important fields

![image](https://github.com/user-attachments/assets/81a8749e-906b-4a8e-bbb6-a32a1e024fc3)

data meta file (one for each session):
* You can manually make this yourself and save as csv
* blockname, name (channel names), and fs (sampling rate) is most important to get right and consistent
* Each row is a specific channel in the data
* 
![image](https://github.com/user-attachments/assets/2f5fe7e3-7760-4e00-9331-7a5304f920fb)


log master file (log_data_fp_tdt.csv):
* This is the master list of sessions to analyze; sessions outlined here must match up with the data in the extracted folder
* copy/paste a new line and edit to register other data to preprocess
* Important parameters:
*   subject: needs to be a substring (streams_data file needs to contain this name) of the streams_data file name
*   blockname: needs to be a substring (streams_data file needs to contain this name) of the streams_data file name
*   signal_id01/control_id: names (signal and isosbestic, respectively) for each channel. Needs to match up with channel names in the streams_data file
*   trim_time_start: amount of time to cut out at the beginning of the session
*   fit_model: NA is default; `huber` if the fit looks bad



![image](https://github.com/user-attachments/assets/78951cf8-a32c-464f-8e37-53a4b6a2e5b4)



