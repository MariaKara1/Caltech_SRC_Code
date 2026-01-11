import csv
import math
import matplotlib.pyplot as plt
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
import astropy.units as u
from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent
csv_path = BASE_DIR / "SRC Obliquity Sheet - Sheet1.csv"

trues = []
obls = []
masses = []
massjs = []
pobls = []
orbs = []
teffs = []
alrobls = []
alrmassjs = []
true_obliquities = []
stellar_temps = []

with open(csv_path, "r") as file:
    reader = csv.reader(file)
    for row in reader:
        if row[13] != "pl_trueobliq":
            true_obliquities.append(float(row[13]))
        if row[2] != "Stellar Temp":
            stellar_temps.append(float(row[2]))

filtered = NasaExoplanetArchive.query_criteria(
    table = "pscomppars", 
    select = "pl_name, pl_projobliq, pl_trueobliq, st_rad, st_vsin, st_rotp, pl_masse, pl_massj, pl_orbper, st_teff, hostname",
    where = "pl_projobliq is not null and pl_trueobliq is null and st_rad is not null and st_vsin is not null and st_rotp is not null")

for row in filtered:
    if row["hostname"] == "TRAPPIST-1":
        row["st_rotp"] = 3.2 * u.day

alrknownts = NasaExoplanetArchive.query_criteria(
    table = "pscomppars",
    select = "pl_name, pl_trueobliq, pl_massj",
    where = "pl_trueobliq is not null")

for row in filtered:
    vr = (2 * math.pi * row["st_rad"] * 696340)/ (row["st_rotp"] * 24* 60 * 60)
    l = ((row["st_vsin"]) / (vr)).value
    if l <= 1:
        trues.append(row)
        inclination = (math.asin(l))
        b = math.acos((math.sin(inclination) * math.cos((row["pl_projobliq"]).to(u.rad).value)))
        true = b * 180/math.pi
        obls.append(true)
        print("The predicted true obliquity for ", row["pl_name"], " is ", true, " degrees. The projected obliquity is ", row["pl_projobliq"], " degrees. The stellar rotation period is ", row["st_rotp"], " days. The stellar radius is ", row["st_rad"], " solar radii. The stellar vsini is ", row["st_vsin"], " km/s. The planet mass is ", row["pl_masse"], " Earth masses or ", row["pl_massj"], " Jupiter masses. The orbital period is ", row["pl_orbper"], " days. The stellar effective temperature is ", row["st_teff"], " K.")
    else:
        trues.append(row)
        l=1
        inclination = (math.asin(l))
        b = math.acos((math.sin(inclination) * math.cos((row["pl_projobliq"]).to(u.rad).value)))
        true = b * 180/math.pi
        obls.append(true)
        print("The predicted true obliquity for ", row["pl_name"], " is ", true, " degrees. The projected obliquity is ", row["pl_projobliq"], " degrees. The stellar rotation period is ", row["st_rotp"], " days. The stellar radius is ", row["st_rad"], " solar radii. The stellar vsini is ", row["st_vsin"], " km/s. The planet mass is ", row["pl_masse"], " Earth masses or ", row["pl_massj"], " Jupiter masses. The orbital period is ", row["pl_orbper"], " days. The stellar effective temperature is ", row["st_teff"], " K.")     

for row in trues:
    masses.append(row["pl_masse"])

for row in trues:
    massjs.append(row['pl_massj'])

for row in trues:
    pobls.append(abs(row["pl_projobliq"].value))

for row in trues:
    orbs.append(row['pl_orbper'].value)

for row in trues:
    teffs.append(row['st_teff'].value)

for row in alrknownts:
    alrobls.append(row['pl_trueobliq'].value)

for row in alrknownts:
    alrmassjs.append(row['pl_massj'])

print("The amount of true obliquities before is ", len(alrobls), ". The amount of true obliquities after is ", len(obls) + len(alrobls), ".")

plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(massjs, obls, alpha=0.6, edgecolors='black')
plt.xlabel("Planet Mass (Jupiter Masses)")
plt.ylabel("True Obliquity (Degrees)")
plt.title("True Obliquity vs Planet Mass")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(massjs, pobls, alpha=0.6, edgecolors='black')
plt.xlabel("Planet Mass (Jupiter Masses)")
plt.ylabel("Projected Obliquity (Degrees)")
plt.title("Projected Obliquity vs Planet Mass")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(orbs, obls, alpha=0.6, edgecolors='black')
plt.xlabel("Orbital Period (Days)")
plt.ylabel("True Obliquity (Degrees)")
plt.title("True Obliquity vs Orbital Period")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(orbs, pobls, alpha=0.6, edgecolors='black')
plt.xlabel("Orbital Period (Days)")
plt.ylabel("Projected Obliquity (Degrees)")
plt.title("Projected Obliquity vs Orbital Period")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(teffs, obls, alpha=0.6, edgecolors='black')
plt.xlabel("Stellar Effective Temperature (K)")
plt.ylabel("True Obliquity (Degrees)")
plt.title("True Obliquity vs Stellar Effective Temperature")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(teffs, pobls, alpha=0.6, edgecolors='black')
plt.xlabel("Stellar Effective Temperature (K)")
plt.ylabel("Projected Obliquity (Degrees)")
plt.title("Projected Obliquity vs Stellar Effective Temperature")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(massjs, obls, alpha=0.6, edgecolors='black')
plt.scatter(alrmassjs, alrobls, color='red', alpha=0.6, edgecolors='black', label='Known True Obliquities')
plt.xlabel("Planet Mass (Jupiter Masses)")
plt.ylabel("True Obliquity (Degrees)")
plt.title("True Obliquity vs Planet Mass")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(alrmassjs, alrobls, color='red', alpha=0.6, edgecolors='black', label='Known True Obliquities')
plt.xlabel("Planet Mass (Jupiter Masses)")
plt.ylabel("True Obliquity (Degrees)")
plt.title("True Obliquity vs Planet Mass")

plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.ylim(-6, 180)
plt.scatter(teffs, obls, alpha=0.6, edgecolors='black')
plt.scatter(stellar_temps, true_obliquities, color = 'red', alpha=0.6, edgecolors='black')
plt.xlabel("Stellar Effective Temperature (K)")
plt.ylabel("True Obliquity (Degrees)")
plt.title("True Obliquity vs Stellar Effective Temperature")

plt.grid(True)
plt.tight_layout()
plt.show()

for row in trues:
    print(row["pl_name"])