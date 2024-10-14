from ast import excepthandler
import csv
from multiprocessing import Value
from pathlib import Path
import easygui

dir_path = easygui.diropenbox()
files_in_dir = [i for i in Path(dir_path).iterdir()]

correct_fpaths = False
while not correct_fpaths:
    fpaths = easygui.multchoicebox('Choose the input files.',
                                   'Files in {}'.format(files_in_dir),
                                   files_in_dir)
    print(f"Got {fpaths} as filepaths to change.")
    try:
        if all([Path(fp).suffix=='.feature' for fp in fpaths]):
            correct_fpaths=True
        else: 
            raise ValueError("Invalid filetype chosen.")
        
    except ValueError:
        easygui.exceptionbox("Please make sure only '.feature' files are chosen.", "ValueError")

# Full/absolute file paths for each file
#fpaths = [r'C:\Users\student\Downloads\05-13-24_NeoV_HIV_TD_B-HCD_ms1.feature']

# Iterate through provided file paths
for p in fpaths:
    try:
        fpath = Path(p)
        print(fpath)
        with open(fpath) as f:
            data = list(csv.reader(f, delimiter='\t'))
    
        # Add more key value pairs to this dictionary as needed. The keys(left of colons) are the old 
        # column names and the values(right of colons) are the new column names.
        from_to = {'Feature_ID': 'ID',
                   'Monoisotopic_mass': 'Mass',
                   'Min_time': 'Time_begin',
                   'Max_time': 'Time_end',
                   'Apex_time': 'Time_apex',
                   'Min_fraction_ID': 'Minimum_fraction_id',
                   'Max_fraction_ID': 'Maximum_fraction_id',
                   'Min_charge': 'Minimum_charge_state',
                   'Max_charge': 'Maximum_charge_state'}

        # The column name change happens here
        for ind, i in enumerate(data[0]):
            if i in list(from_to.keys()):
                data[0][ind] = from_to[i]
        
        # Create a new file name for the new updated file with an assignable modifier. 
        fname_modifier = '_updated'
        new_fpath = fpath.parent / (fpath.stem + fname_modifier + fpath.suffix)

        # Save to same folder as the input file.
        with open(new_fpath, 'w', newline='') as output:
            writer = csv.writer(output, delimiter='\t')
            for line in data:
                writer.writerow(line)
    except:
        print('Something went wrong with file path {}'.format(p))


