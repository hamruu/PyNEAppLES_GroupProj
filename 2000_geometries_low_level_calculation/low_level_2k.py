import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("C://Users//marij//OneDrive//Desktop//sci-comp group//2000_low.tsv", sep="\t")

x = df.iloc[:, 0] 
y = df.iloc[:, 1] 

#plot of whole spectrum
plt.figure(figsize=(10, 6))
plt.plot(x, y, marker='o', linestyle='-', markersize=1)
plt.xlabel("Energy (eV)")
plt.ylabel("Cross-section (cm$^2$)")
plt.title("Spectrum of 2000 Geometries Low-Level Calculation")
plt.grid(True)
plt.show()


#x-range (showing only the first excited state)
x_min, x_max = 3.0, 5.5

#filter rows where x is between 3 and 5.5
df_filtered = df[(df.iloc[:, 0] >= x_min) & (df.iloc[:, 0] <= x_max)]

#check if there's data in the range
if df_filtered.empty:
    print(f"No data in the range [{x_min}, {x_max}]")
else:
    x = df_filtered.iloc[:, 0]
    y = df_filtered.iloc[:, 1]

    # Plot the filtered data
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='-', markersize=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("Cross-section (cm$^2$)")
    plt.title("Spectrum of 2000 Geometries Low-Level Calculation for 1st Excited State")
    plt.grid(True)
    plt.show()