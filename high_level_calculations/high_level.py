import pandas as pd
import matplotlib.pyplot as plt

#list of file paths and corresponding labels
files = [
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//high level//spectrum-1-geom.tsv", "Spectrum of 1 Representative Geometry"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//high level//spectrum-3-geoms.tsv", "Spectrum of 3 Representative Geometries"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//high level//spectrum-5-geoms.tsv", "Spectrum of 5 representative Geometries"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//high level//spectrum-10-geoms.tsv", "Spectrum of 10 representative Geometries"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//high level//spectrum-30-geoms.tsv", "Spectrum of 30 representative Geometries"),
    ("C://Users//marij//OneDrive//Desktop//sci-comp group//high level//spectrum-30-geoms.tsv", "Spectrum of 50 representative Geometries")    
]

#plot each file in its own figure
for file_path, label in files:
    df = pd.read_csv(file_path, sep="\t")
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='-', markersize=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Cross-section (cm$^2$)")
    plt.title(label)
    plt.grid(True)
    plt.show()

#plot all files in a single graph
plt.figure(figsize=(10, 6))
for file_path, label in files:
    df = pd.read_csv(file_path, sep="\t")
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]
    plt.plot(x, y, marker='o', linestyle='-', markersize=1, label=label)   

plt.xlabel("Energy (eV)")
plt.ylabel("Cross-section (cm$^2$)")
plt.title("Combined Spectrum for High Level Calculations")
plt.legend()
plt.grid(True)
plt.show()