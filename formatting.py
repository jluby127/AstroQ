import numpy as np
import pandas as pd
import star

def theTTP(filename):

    targets = pd.read_csv(filename)

    if list(targets.columns) != ['Starname', 'RA', 'Dec', 'Exposure Time', 'Exposures Per Visit', 'Visits In Night', 'Intra_Night_Cadence', 'Priority']:
        print(list(targets.columns))
        print("Error: column names not correct.")
        print("Column names must be this format and this order: ['Starname', 'RA', 'Dec', 'Exposure Time', 'Exposures Per Visit', 'Visits In Night', 'Intra_Night_Cadence', 'Priority'].")
        return
    else:
        print("Building Star objects:")
        all_star_objects = []
        for t, row in targets.iterrows():
            starobj = star.star(row, t)
            #starobj.printStar()
            all_star_objects.append(starobj)
        return all_star_objects
