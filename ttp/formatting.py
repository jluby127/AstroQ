import numpy as np
import pandas as pd
import ttp.star #for testing, use just star

def theTTP(filename):
    """Read in a .csv file of targets and convert them to star objects

    Args:
        filename (string): Path to .csv file

    Returns:
        all_star_objects (list): List of star objects, one for each row
            in the .csv file
    """

    targets = pd.read_csv(filename)

    required_columns = ['Starname', 'RA', 'Dec', 'Exposure Time', 'Exposures Per Visit', 'Visits In Night', 'Intra_Night_Cadence', 'Priority']
    optional_columns = ['First Available', 'Last Available']
    
    # Check if the first 8 columns match the required columns
    if list(targets.columns)[:8] != required_columns:
        print(list(targets.columns))
        print("Error: column names not correct.")
        print("Column names must be this format and this order: ['Starname', 'RA', 'Dec', 'Exposure Time', 'Exposures Per Visit', 'Visits In Night', 'Intra_Night_Cadence', 'Priority'] or can include ['First Available', 'Last Available'] at the end.")
        return
    else:
        print("Building Star objects:")
        all_star_objects = []
        for t, row in targets.iterrows():
            starobj = star.star(row, t)
            # starobj.printStar()
            all_star_objects.append(starobj)
        return all_star_objects


def readInputs(filename):
    '''
    Read in the inputs from the special formatted file

    filename (str) - the path and filename which holds the input information

    Returns:
        result_dict (dictionary) - the processed input information

    '''
    result_dict = {}
    # Open the file in read mode
    with open(filename, 'r') as file:
        for line in file:
            # Strip any extra whitespace (like newlines) and split by colon
            line = line.strip()
            if ':' in line:
                key, value = line.split(' : ', 1)  # Split only at the first colon
                # important that there be one space after the key name and one space after the colon
                result_dict[key] = value

    return result_dict
