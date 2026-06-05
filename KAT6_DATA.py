# Script written by : sbernard@kat6.org
# TODO : Divide the KAT6A and KAT6B patients, applying some filtering to the data, for graphs, and other conclusions.

import pandas as pd
import matplotlib.pyplot as plt

def kat6_var_diagnosis(data):
    print("Number of Unique Participants, based on ID", data["Participant ID"].nunique())

    # Filtering by genetic variant
    kat6_variants = {}
    for dx, dx_data in data.groupby("dx"):
        kat6_variants[dx] = dx_data["Participant ID"].unique().tolist()
        print(f"Type of KAT6 variant : " f'{dx}: {dx_data["Participant ID"].nunique()} participants')

    # Filtering by diagnosis type
    diagnosis_variants = {}
    for group_variant, group_data in data.groupby("genetic_dx"):
        diagnosis_variants[group_variant] = group_data['Participant ID'].unique().tolist()
        print(f"Type of diagnosis : '{group_variant}': {group_data['Participant ID'].nunique()} participants")

    # Graph
    variant_diagnosis=data.groupby(["dx","genetic_dx"])["Participant ID"].nunique().unstack(fill_value=0)
    ax=variant_diagnosis.plot(kind="bar", figsize=(15,10), colormap="viridis")
    for container in ax.containers:
        ax.bar_label(container, label_type="edge",padding=0.1)

    plt.title("Participants by KAT6 Variant and Diagnosis Type")
    plt.xlabel("KAT6 Variant (dx)")
    plt.ylabel("Number of Participants")
    plt.xticks(rotation=0)
    plt.legend(title="Diagnosis Type (genetic_dx)")
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    kat6_var_diagnosis(data=pd.read_excel(r"C:\Users\sophi\Documents\KAT6\All Questions downloaded May 18.xlsx"))

#TODO : verbal-sex graph ? where can i find the verbal data ?