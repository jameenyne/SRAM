#!/usr/bin/env python3

# Import required libraries
import getopt, numpy as np, pylab as plt, pandas as pnd, os, sys, scipy, pyDOE
import seaborn as sns
from sklearn import linear_model
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pickle

#Loading the samples
pickle_file_path = "sample.out"
# Read the pickle file
with open(pickle_file_path, "rb") as f:
    PARAMS = pickle.load(f).to_numpy()

# Function to read waveform data from a SPICE file
def read_spice_file(wavefile) :
    f = open(wavefile, "r")
    lines = f.readlines()
    f.close()
    READING = False
    for line in lines :
        n = line.split()
        if len(n) == 0 or n[0] == "End" :
            continue
        if n[0] == "Index" :
            DATA = np.zeros((len(lines)-2, len(n)-1))
            continue
        idx = int(n[0])
        DATA[idx,:] = list(float(x) for x in n[1:])
    return DATA

# Function to write circuit parameters to a SPICE file
def write_spice_file(cktfile, d0, d1, d2, d3, d4) :
    f = open(cktfile, "w")
    f.write(f"""* sram 32nm cross-section
.include sram.lib
.param vdd=0.7
* make the DC voltage sources
VDD	vdd	0	{vdd}
VWL1	vwl1	0	0.0	pulse(0.0 {vdd} 0.2n 0.05n 0.05n 0.2n 1.2n)
VWL2	vwl2	0	0.0	pulse(0.0 {vdd} 0.8n 0.05n 0.05n 0.2n 1.2n)
Vpchg	pchg	0	0.0	pulse({vdd} 0.0 0.05n 0.05n 0.05n 0.05n 0.6n)
Venable	enable	0	0.0	pulse(0.0 {vdd} 0.3n 0.025n 0.025n 0.05n 0.6n)
* the cross-section
x0 vwl1 vwl2 pchg enable out_t out_c vdd xsection PARAMS: D0={d0} D1={d1} D2={d2} D3={d3} D4={d4}
* simulation options
.tran 0.01n 1.2n
.print tran v(pchg) v(enable) v(out_t) v(out_c)
.end""")
    f.close()

# Function to plot waveforms from a SPICE file
def plot_waveform(wavefile):
    DATA = read_spice_file(wavefile)
    axs[0].plot(DATA[:, 0], DATA[:, 1], lw=0.5, color = "blue", label="PCHG")
    axs[1].plot(DATA[:, 0], DATA[:, 2], lw=0.5, color = "brown", label="EN")
    axs[2].plot(DATA[:, 0], DATA[:, 3], lw=0.5, color = "red", label="out_t")
    axs[3].plot(DATA[:, 0], DATA[:, 4], lw=0.5, color = "purple", label="out_c")
    axs[0].set_ylabel("V(PCHG)")
    axs[1].set_ylabel("V(Sense)")
    axs[2].set_ylabel("V(BLT)")
    axs[3].set_ylabel("V(BLC)")
    axs[3].set_xlabel("Time")

fig, axs = plt.subplots(4, sharey=True)
fig.suptitle("Waveforms from Xyce")
plt.tight_layout()

# Function to calculate the crossing point of a voltage waveform
def get_cross(time, volt, vdd, tmax) :
    where = time < tmax
    time = time[where]
    volt = volt[where]
    where = np.abs(volt-0.5*vdd) < 0.4*vdd
    if np.sum(where) < 2 :
        return 0.0
    b,a = np.polyfit(time[where], volt[where], 1)
    crossing = (vdd*0.5 - a)/b
    return crossing

# Main simulation function
def run(P, vdd) :
    write_spice_file(crkt_file, P[0], P[1], P[2], P[3], P[4])
    res = os.system(f"C:\\\"Program Files\"\\\"Xyce 7.6 NORAD\"\\bin\\xyce  {crkt_file} > nul 2>&1")
    DATA = read_spice_file(wave_file)
    DATA[:,0] *= 1e12
    plot_waveform(wave_file)
    D0 = get_cross(DATA[:,0],DATA[:,3], vdd, 1200)
    D1 = get_cross(DATA[:,0],DATA[:,4], vdd, 600)
    D0 = D0-912.6 if D0 > 912.6 else None
    D1 = D1-312.56 if D1 > 312.56 else None
    return D0,D1

# Function to display usage information
def usage() :
    print("""usage:
sample.py [flags]
    -h        give this help
    -n <N>    set the number of samples [1000]
    -s <max>  set the range of LHS sampling in sigma [3.5]
    -l        do latin hypercube sampling [simple MC]
""")
    sys.exit()

# Process command-line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:], "hn:s")
except getopt.GetoptError as err:
    raise Exception(str(err))

# Initialize variables and set default values
NSAM=1000
SMAX=4
crkt_file = "x.ckt"
wave_file = "x.ckt.prn"

# Process command-line flags
for o, a in opts :
    if o == "-n" :
        NSAM = int(a)
    elif o == "-s" :
        SMAX = float(a)
    elif o == "-h" :
        usage()

# Set initial values for variables
vdd = 0.7
NVAR = 5
DATA = np.zeros((NSAM,7))

# Loop over samples and perform simulations
for i in range(NSAM) :
    DATA[i,0:NVAR] = PARAMS[i,:]
    DATA[i,NVAR],DATA[i,NVAR+1] = run(PARAMS[i,:], vdd)

# Create a DataFrame from the simulation results
DF = pnd.DataFrame(DATA, columns=("D0", "D1", "D2", "D3", "D4", "Delay0", "Delay1"))

# Create a pair plot to visualize the relationships between input parameters and delay values
DF2 = DF.iloc[:, :5]
DF3 = DF.iloc[:, 5:7]
pairplot = sns.pairplot(DF,x_vars=DF2,y_vars=DF3)
pairplot.fig.suptitle("Latin hypercube sampling (LHS)")
pairplot.fig.tight_layout()

# Perform linear regression analysis for Delay 0
X = DF[["D0","D1","D2","D3","D4"]]
Y0 = DF["Delay0"]
Y1 = DF["Delay1"]
fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(15,7))
regr = linear_model.LinearRegression()

regr.fit(X, Y0)

Y0_HAT = regr.predict(X)
L = axes[0].plot(Y0, Y0_HAT, '.')
ymin,ymax = np.amin(Y0)*0.9,np.amax(Y0)*1.1
axes[0].plot((ymin,ymax), (ymin,ymax), '--', color=L[-1].get_color(), lw=0.5)
axes[0].set_xlabel("Delay 0")
axes[0].set_ylabel("Linear Model")

#Sensitvity analysis for Delay 0
BX = inset_axes(axes[0], width="40%", height="40%", bbox_to_anchor=(0, -0.53, 1, 1), bbox_transform=axes[0].transAxes)
coefs = regr.coef_ / np.linalg.norm(regr.coef_)
BX.bar(range(5), coefs)
BX.set_xticks(range(5))

# Perform linear regression analysis for Delay 1
regr.fit(X, Y1)
Y1_HAT = regr.predict(X)
COEF = regr.coef_
L = axes[1].plot(Y1, Y1_HAT, '.')
ymin,ymax = np.amin(Y1)*0.9,np.amax(Y1)*1.1
axes[1].plot((ymin,ymax), (ymin,ymax), '--', color=L[-1].get_color(), lw=0.5)
axes[1].set_xlabel("Delay 1")
axes[1].set_ylabel("Linear Model")

#Sensitvity analysis for Delay 1
BX = inset_axes(axes[1], width="40%", height="40%", bbox_to_anchor=(0, -0.53, 1, 1), bbox_transform=axes[1].transAxes)
coefs = regr.coef_ / np.linalg.norm(regr.coef_)
BX.bar(range(5), coefs)
BX.set_xticks(range(5))
fig.suptitle("Sensitivity of Delay 0 and Delay 1")


# Plot histograms of Delay 0 and Delay 1
fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

mean0 = np.mean(DF["Delay0"])
stdv0 = np.std(DF["Delay0"])
n, bins, patches = axes[0].hist(DF["Delay0"], 30, density=True, facecolor='blue', alpha=0.5)
y = scipy.stats.norm(mean0,stdv0).pdf(bins)
l = axes[0].plot(bins, y, 'r--', linewidth=1)
axes[0].set_title(f"Delay 0 {mean0:.4e} {stdv0:.4e}")

mean1 = np.mean(DF["Delay1"])
stdv1 = np.std(DF["Delay1"])
n, bins, patches = axes[1].hist(DF["Delay1"], 30, density=True, facecolor='blue', alpha=0.5)
y = scipy.stats.norm(mean1,stdv1).pdf(bins)
l = axes[1].plot(bins, y, 'r--', linewidth=1)
axes[1].set_title(f"Delay 1 {mean1:.4e} {stdv1:.4e}")

# Calculate and plot the failure probability from Delay>40 till Delay=60
D = np.linspace(41, 60, 200)
p0 = 1.0 - scipy.stats.norm.cdf(D, mean0, stdv0)
p1 = 1.0 - scipy.stats.norm.cdf(D, mean1, stdv1)
FAIL = p0 + p1 - p0*p1
fig,axs = plt.subplots(figsize=(10,8))
axs.plot(D, FAIL, '.-')
axs.set_yscale("log")
axs.set_xlabel("Max Delay (ps)")
axs.set_ylabel("Failure Probability")
axs.set_title("Max Delay vs Failure Probability for Sense Amplifier")
fail_len = len(FAIL)
print('Fail rate: ', format(FAIL[100],".6f"))

# Plot the relationship between SA Fail Probability and Array Fail Probability
fig,ax = plt.subplots(figsize=(10,8))
for n in FAIL:
    p_ARR_succ = (1-FAIL)**64 #Calculates array success probability for 64 sense amp in a 8K array 
    p_ARR_fail = 1-p_ARR_succ #Array fail probability
    ax.plot(FAIL, p_ARR_fail, '-')

x_marker = FAIL[100]
y_marker = p_ARR_fail[100]
y_intercept = 1-(1-FAIL[100])**64
ax.text(x_marker, p_ARR_fail[199], f'SA Fail Probability = {round(x_marker, 6)}', ha='left')
ax.text(FAIL[199], y_marker, f'Array Fail Probability = {round(y_marker, 6)}', va='bottom')
ax.plot([x_marker, x_marker], [0, y_intercept], color='red', linestyle='--')
ax.plot([0, x_marker], [y_marker, y_marker], color='blue', linestyle='--')
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel("SA Fail Probability")
ax.set_ylabel("Array Fail Probability")
ax.set_title("SA Fail Probability vs Array Fail Probability")
plt.tight_layout()
plt.show()