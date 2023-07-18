import numpy as np
import matplotlib.pyplot as plt
import time
import subprocess
import pandas as pd
import random
from scipy.optimize import curve_fit

# choose R (reasonable) or G (gubser)
LAMBDA = "G"

FILENAME = "out/output_" + LAMBDA + "b.csv"

# length of csv file (N) obtained using bash command
N = int(subprocess.check_output("wc -l " + FILENAME, shell = True).split()[0])

# load csv file into pandas dataframe
df = pd.read_csv(FILENAME)
print(f"Filename: {FILENAME}")
print("Excerpt of collected data:\n",df)

# cuts based on transverse momenta of leading and secondary quark
# cut1 - Leading Quark must be between 10 and 20 GeV

LeadCuts = input("Please enter the Lead Quark Lower and Upper Bounds:\n").strip()
if LeadCuts == "":
    # default values
    leadLower = 10
    leadUpper = 20
else:
    leadLower, leadUpper = LeadCuts.split("-")
    leadLower = int(leadLower)
    leadUpper = int(leadUpper)

SecondCuts = input("Please enter the Secondary Quark Lower and Upper Bounds:\n").strip()
if SecondCuts == "":
    # default values
    secondLower = 5
    secondUpper = 10
else:
    secondLower, secondUpper = SecondCuts.split("-")
    secondLower = int(secondLower)
    secondUpper = int(secondUpper)

# determine runtime
start_time = time.time()

# given that the angular separation is between 0 and 2*pi
# cuts1 - Leading Quark must be between 10 and 20 GeV
cuts1 = df[(df["quarkLead_PT"] > leadLower) & (df["quarkLead_PT"] < leadUpper )]

# cuts2 - Secondary Quark must be between 5 and 10 GeV
cuts2 = cuts1[(cuts1["quarkSecondary_PT"] > secondLower) & (cuts1["quarkSecondary_PT"] < secondUpper)]

# setting final dataframe 
finalData = cuts2
#finalData = df

# binning procedure
nBins = input("How many bins do you want? (preferably odd)\n").strip()
if nBins == "":
    # default value
    nBins = 175
else:
    nBins = int(nBins)
minAngle = min(finalData['deltaAngle'])
maxAngle = max(finalData['deltaAngle'])

countsData, edges = np.histogram(finalData['deltaAngle'], bins=nBins, range = (0,2*np.pi))

# create an array at centre of each bin
centres = (edges[1:] + edges[:-1])/2

# getting mean and FWHM
#print(countsData, centres)
mean_index = np.argmax(countsData)
mean = centres[mean_index]

half_max = max(countsData)/2
# need to find element in countsData closest to half_max
closest = 0
difference = float('inf')
for i in range(0,len(countsData)):
    if (np.abs(countsData[i]-half_max)<difference):
        difference = np.abs(countsData[i]-half_max)
        closest = i
angle_fwhm = centres[closest]
fwhm = np.abs(centres[closest] - mean)*2

# finds closest bin to pi for comparison to mean centre
# (since bins aren't always going to be exactly at pi)
closest = 0
difference = float('inf')
for i in range(0,len(centres)):
    if (np.abs(centres[i]-np.pi)<difference):
        difference = np.abs(centres[i]-np.pi)
        closest = i
nearestBin = centres[closest]


# priting info to console
print(f"Size of Dataset Post-Cuts: {finalData.size:,d}")
print("+--------------------------------+")
print("|    Distribution Information    |")
print("+--------------------------------+")
print(f"| Mean: {mean:.6f} rad             |")
print(f"| FWHM: {fwhm:.6f} rad             |")
print("+--------------------------------+")
print("|     Additional Information     |")
print("+--------------------------------+")
print(f"| Closest to π: {nearestBin:.6f} rad     |")


# Uncertainty calculation (using Poisson approximation ~ sqrt(N) for now)
unc_countsData = np.zeros(len(countsData))
for i in range(0,len(countsData)):
    if countsData[i] != 0:
        unc_countsData[i] = 1/np.sqrt(countsData[i])
    else:
        unc_countsData[i] = 0 

# plotting histogram
plt.figure(0)
plt.title(f"Leading Quark Between {leadLower} and {leadUpper} GeV\nSecondary Quark Between {secondLower} and {secondUpper} GeV")
#plt.hist(finalData['deltaAngle'], bins = nBins, range = [minAngle,maxAngle], color = "green", alpha=0.6)
plt.errorbar(centres, countsData, yerr = unc_countsData, ecolor = "r", color = "k", fmt = ".", label = "Binned Distribution")
plt.xlabel("Angular Separation (rad)")
plt.ylabel("Number of Pairs (counts)")
plt.xlim(0,2*np.pi)
plt.xticks([0, np.pi/2, np.pi, (3/2)*np.pi, 2*np.pi], ['0', 'π/2', 'π', '3π/2', '2π'])

# Alternative Uncertainty Calculations

# split dataframe into 100 equal parts
numDivisons = 10
dfArray = np.array_split(finalData, numDivisons)

# making each df into a histogram
countsArray = []
edgesArray = []
centresArray = []

# testing
x = random.randint(0,9)

for i in range(0, len(dfArray)):
    c, e = np.histogram(dfArray[i]['deltaAngle'], bins=nBins, range = (0,2*np.pi))
    countsArray.append(c)
    edgesArray.append(e)
    # create an array at centre of each bin
    centre = (edgesArray[i][1:] + edgesArray[i][:-1])/2
    centresArray.append(centre)

meanArray = []
uncertaintyArray = []
# calculating mean value for each bin
lenArray = len(countsArray)
for j in range(0, nBins):
    sumCounts = 0
    for k in range(0, lenArray):
        sumCounts += countsArray[k][j]
    m = sumCounts/lenArray
    meanArray.append(m)

# calculating uncertainty for each bin
for j in range(0, nBins):
    varCounts = 0
    for k in range(0, lenArray):
        varCounts += (countsArray[k][j] - meanArray[j])**2
    var = varCounts/(lenArray-1)
    std = np.sqrt(var)
    uncertaintyArray.append(std/np.sqrt(lenArray))

plt.figure(2) 
plt.title(f"Plot when dataset is split into {numDivisons} subsets")
plt.scatter(centresArray[x], countsArray[x], marker=".", color = "green", label = f"Example of Subset #{x}")
plt.errorbar(centresArray[x], meanArray, yerr = uncertaintyArray, ecolor = "red", color = "k", fmt = ".", label = "Mean Values")
plt.xlabel("Angular Separation (rad)")
plt.ylabel("Number of Pairs (counts)")
plt.xlim(0,2*np.pi)
plt.xticks([0, np.pi/2, np.pi, (3/2)*np.pi, 2*np.pi], ['0', 'π/2', 'π', '3π/2', '2π'])
plt.legend()

# fitting procedure for gaussian to find mean and fwhm
# need fitting dataset that doesn't include 0 count bins
MAXIMUM_COUNTS = max(countsData)
THRESHOLD = int(MAXIMUM_COUNTS*0.35)
fit_countsData = []
fit_centres = []
for i in range(0,len(countsData)):
    if countsData[i] >= THRESHOLD:
        fit_countsData.append(countsData[i])
        fit_centres.append(centres[i])

print(f"| Maximum: {MAXIMUM_COUNTS:5,d} counts          |")
print(f"| Threshold: {THRESHOLD:5,d} counts        |")


plt.figure(0)

# gaussian function 
def g(x, A, mu, sigma):
    return A*np.exp(-1*((x-mu)**2)/(2*sigma**2))

# estimates
A0 = 6000
mu0 = 3.1415926
sigma0 = 0.2

p0 = [A0, mu0, sigma0]
name = ["A0", "mu0", "sigma0"]

xmodel = np.linspace(minAngle,maxAngle,10000)
ystart = g(xmodel,*p0)

# perform the fit
udata = np.ones(len(fit_countsData))
popt, pcov = curve_fit(g,fit_centres,fit_countsData,p0,sigma=udata,absolute_sigma=True)

# plotting final fitted line
yfit = g(xmodel,*popt)
plt.plot(xmodel,yfit,'-k', label = 'Fitted Gaussian')

# printing fitted parameters
fitted_fwhm = popt[2]*2*np.sqrt(2*np.log(2))
print("+--------------------------------+")
print("|      Gaussian Information      |")
print("+--------------------------------+")
print(f"| Mean: {popt[1]:.6f} rad             |")
print(f"| FWHM: {fitted_fwhm:.6f} rad             |")
print("+--------------------------------+")

# try fit Breit-Wigner function to distribution
def f(x, A, x0, P):
    return P * (A / (((x-x0)**2 + A**2) * np.pi))

# estimates
A0 = 0.2
x0 = 3.1415926
P0 = 5000

p0 = [A0, x0, P0]
# doing fit
#norm_factor = np.trapz(fit_countsData, fit_centres)
popt2, pcov2 = curve_fit(f, fit_centres, fit_countsData, p0)

# plot fitted curve
y_val = f(xmodel, *popt2)
plt.plot(xmodel, y_val, label = "Fitted Cauchy", color = "g" )

# printing fitted parameters
fitted_fwhm2 = popt2[0]*2
print("|       Cauchy Information       |")
print("+--------------------------------+")
print(f"| Mean: {popt2[1]:.6f} rad             |")
print(f"| FWHM: {fitted_fwhm2:.6f} rad             |")
print("+--------------------------------+")

# print runtime measured to console (to 3 decimal places)
print(f"Runtime of {time.time() - start_time:.3f} seconds")

# display all plotted graphs
plt.legend()
plt.show()