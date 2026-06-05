#Script to divide the patients in age groups and gender, produce a graph for that data
# output a "master file" type of excel that gives us space to keep track of added information patients send back

import pandas as pd
import matplotlib.pyplot as plt


if __name__ == '__main__':
    fonction(data=pd.read_excel(r"C:\Users\sophi\Documents\KAT6\All Questions downloaded May 18.xlsx"))