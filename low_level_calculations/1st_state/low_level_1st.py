import pandas as pd
import matplotlib.pyplot as plt

#list of file paths and corresponding labels
files = [
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//low level//spectrum-1-geom-low.tsv", "Spectrum of 1 Representative Geometry for 1st Excited State"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//low level//spectrum-3-geom-low.tsv", "Spectrum of 3 representative Geometries for 1st Excited State"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//low level//spectrum-5-geom-low.tsv", "Spectrum of 5 representative Geometries for 1st Excited State"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//low level//spectrum-10-geom-low.tsv", "Spectrum of 10 representative Geometries for 1st Excited State"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//low level//spectrum-30-geom-low.tsv", "Spectrum of 30 representative Geometries for 1st Excited State"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//low level//spectrum-50-geom-low.tsv", "Spectrum of 50 representative Geometries for 1st Excited State"),
]

#x-range (showing only the first excited state)
x_min, x_max = 3.0, 5.5

#plot each file in its own figure, showing only x in [3, 5.5]
for file_path, label in files:
    df = pd.read_csv(file_path, sep="\t")
    
    #filter rows where x is between 3 and 5.5
    df_filtered = df[(df.iloc[:, 0] >= x_min) & (df.iloc[:, 0] <= x_max)]
    
    #if there's no data in that range, df_filtered may be empty
    if df_filtered.empty:
        print(f"No data for {label} in the range [{x_min}, {x_max}]")
        continue
    
    x = df_filtered.iloc[:, 0]
    y = df_filtered.iloc[:, 1]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='-', markersize=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Cross-section (cm$^2$)")
    plt.title(label)
    plt.grid(True)
    plt.show()

#plot all files on a single graph
plt.figure(figsize=(10, 6))
for file_path, label in files:
    df = pd.read_csv(file_path, sep="\t")
    df_filtered = df[(df.iloc[:, 0] >= x_min) & (df.iloc[:, 0] <= x_max)]
    
    #skip if there's no data in that range
    if df_filtered.empty:
        print(f"No data for {label} in the range [{x_min}, {x_max}]")
        continue
    
    x = df_filtered.iloc[:, 0]
    y = df_filtered.iloc[:, 1]
    plt.plot(x, y, marker='o', linestyle='-', markersize=1, label=label)

plt.xlabel("Energy (eV)")
plt.ylabel("Cross-section (cm$^2$)")
plt.title(f"Combined Spectrum for the 1st Excited State for Low Level Calculations")
plt.grid(True)
plt.legend()
plt.show()