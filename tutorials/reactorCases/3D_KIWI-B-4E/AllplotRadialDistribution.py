"""
Script to plot radial distributions from multiple files. Plot
- Neutron fluxes (group 0 and 1)
- Power density

Thomas Guilbaud, EPFL, 13/08/2023
"""

#==============================================================================*
# Imports
#==============================================================================*

import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


#==============================================================================*
# Constants
#==============================================================================*

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


#==============================================================================*
# Reactor Parameters
#==============================================================================*

dzCore = 1.32 # m

rCore = 0.43180
roReflector = 0.508
roReflectorSystem = 0.62547

#==============================================================================*
# Usage
#==============================================================================*

if (len(sys.argv) <= 1):
    print(f"\n  Usage: python3 {sys.argv[0]} path/to/data_<>.xy path2/to/data_<>.xy ...\n")
    sys.exit(1)


#==============================================================================*
# Usefull Functions and Classes
#==============================================================================*

def listModifier(l: list, scale: float=1, offset: float=0) -> list:
    return([scale*e+offset for e in l])


def doubleArrow(ax, x, y, dx, head_width, color):
    ax.arrow(
        x=x, y=y, dx=dx, dy=0, length_includes_head=True,
        head_width=head_width, color=color
    )
    ax.arrow(
        x=dx, y=y, dx=-dx, dy=0, length_includes_head=True,
        head_width=head_width, color=color
    )


def plotRefFission(ax):
    """ Place experimental measures of relative fission distribution in K from LA-3327-MS """
    ax.scatter(
        [r/1e2 for r in [1.54490959211875, 1.54859558747497, 1.55072958478647, 1.53300099173709,
         1.55519157916506, 1.54767035787139, 3.21984378965773, 3.71660254843531,
         4.87173678545577, 5.5728220560449, 6.48921720922363, 6.7811420722135,
         7.52389252108029, 8.23739377602722, 8.65683617067279, 9.83163899829854,
         10.550452861937, 11.5250777879093, 11.8499378401752, 13.3057269291872,
         13.7233039415674, 14.8955847723705, 15.5039530828453, 16.5278837928538,
         16.8684578253227, 18.250885160605, 18.7244833331767, 19.6946462647703,
         20.4872098816496, 21.6063795485951, 21.9047064035195, 23.3681360982901,
         23.7472713129471, 24.7966906062287, 25.5999988249561, 26.6094989377603,
         26.963727568411, 28.4020566794198, 28.8538522640747, 29.9242981462507,
         30.6680185938954, 31.7420310869626, 31.9973347653202, 33.099820453285,
         33.4715836772295, 33.5941318305305, 33.789071738783, 33.9253938747309,
         34.2807865039152, 34.5164215916675, 34.5575047706784, 34.4636984273212,
         34.7139581120334, 34.8473702516474, 35.0810653418438, 35.3506056176501,
         35.611311443048, 35.6085954464697, 35.7205929976781, 36.2824461359855,
         36.2343789657733, 36.4328108696265, 36.2856247473655, 36.3592551161413,
         36.4303634880944, 36.5855782156253, 36.5593136333299, 36.8159454638604,
         36.8570286428712, 37.1328515261471, 37.2967962426795, 37.2483410729561,
         37.5268053375196, 37.5566365307062, 37.5591585275289, 38.0017615083804,
         38.0319807010782, 38.2541551904041, 38.2767337773433, 38.2800317731883,
         38.1234590473684, 38.1263690437022, 38.2862397653672, 38.1046858402504,
         38.5770752451142, 39.0858619887384, 39.5102289925643, 39.4488952236814,
         39.4642958196636, 39.4765924195565, 38.9309606454281, 38.6970715554762,
         38.5825520074451, 38.7026975483883, 39.6979162945694, 39.7389994735803,
         39.9006161930456, 40.0141657422988, 39.870665615394, 40.0846921149851,
         40.000779759163, 40.1558004869382, 40.2518154428976, 40.3685137574145,
         40.7865533845966, 40.7910153789752, 41.0244418541253, 41.3038761174668,
         41.3132627210258, 41.3111287237143, 41.4200222788332, 41.3645831179086,
         41.4060542964307, 41.5472711954427, 41.5715957801822, 41.8793092386655,
         41.8769812415984, 41.8540146551481, 41.848388662236, 41.8263920745636,
         42.0410005734214, 42.0564011694037, 41.9345096306602, 42.3674872390227,
         42.3846338328053]],
        [1.02526567743634, 1.01805815370828, 1.01388537681309, 1.01008998486543,
         1.00516047966861, 0.94294425278172, 1.01936264832971, 1.00954868128735,
         0.981590860664855, 0.995317724049546, 0.972647209505635, 0.978746079270031,
         0.98792473921858, 0.977373627940351, 0.964896423792916, 0.936940561577002,
         0.954462816909352, 0.933315551001607, 0.951936668264726, 0.913007266471766,
         0.942639133036908, 0.9196147344244, 0.922330848507067, 0.920156625524461,
         0.908050931104823, 0.897187062296129, 0.89419931722113, 0.881776948457872,
         0.870474984097738, 0.87438004681375, 0.867960585892578, 0.852552430460899,
         0.841968026270849, 0.828415461071577, 0.834565249406994, 0.822146797456891,
         0.821802705421182, 0.817014205497952, 0.818197278911564, 0.801991268640113,
         0.809273211818121, 0.824554658344221, 0.825338804337948, 0.861867395507258,
         0.865698038773316, 0.895300529239793, 0.898734203492544, 0.901403315817344,
         0.898783163656987, 0.899565351244136, 0.896155177870318, 0.92573612586444,
         0.936383199064978, 0.9447424617014, 0.949318082829631, 0.960725801144962,
         0.835562078355063, 0.840872885312577, 0.852644279729395, 0.869392181100906,
         0.886458715222616, 0.88306420910142, 0.901638324606673, 0.911509085439795,
         0.926311309876322, 0.930499754024133, 0.943395469657231, 0.826197957303602,
         0.822787783929785, 0.860372152085154, 0.8474901452981, 0.865315366128026,
         0.859273290154385, 0.877865031318838, 0.872933567715432, 0.853630141600629,
         0.871463196056865, 0.860104437905977, 0.854416246000933, 0.847967408981095,
         0.846434368312041, 0.840744218000419, 0.835828421649636, 0.806220055963426,
         0.882518988710177, 0.849185733713107, 0.865540582884466, 0.87008682791403,
         0.878434340110986, 0.892851345973672, 0.805922965685583, 0.801726687911461,
         0.794886757098048, 0.790725730642325, 0.84469432406772, 0.841284150693902,
         0.832954264156145, 0.841690911740098, 0.853057503517297, 0.85763116623895,
         0.867865603173558, 0.872433390675477, 0.876995302957663, 0.802651839178785,
         0.869840656207208, 0.861115759062722, 0.827755086373563, 0.819816293629381,
         0.839923449803689, 0.844096226698878, 0.861937114781426, 0.893418696359243,
         0.88924983627721, 0.805424747052206, 0.796322464960189, 0.810010160213325,
         0.814562280462622, 0.821009159075882, 0.832010116345017, 0.836560278187737,
         0.839995910847065, 0.848343423044021, 0.855918539686717, 0.893897918448816,
         0.898831340458799
        ],
        label="Experimental data",
        marker='x'
    )



class Figure:
    """
    Figure object to simplify the creation of graphs.
    """
    def __init__(
            self,
            title: str,
            yLabel: str,
            isYlimZero: bool=False,
            isYlog: bool=False
        ) -> None:
        self.fig, self.ax = plt.subplots(figsize=(9, 5))
        # self.fig.subplots_adjust(left=0.15, right=0.85)
        self.ax.set_title(title)
        self.ax.set_xlabel("Reactor radial position [m]")
        self.ax.set_ylabel(yLabel)
        self.ax.grid()
        self.isYlimZero = isYlimZero

        if (isYlog):
            self.ax.set_yscale('log')

        self.allYList = []

    def plot(self, x: list, y: list, label: str, ls: str, color: str=None) -> None:
        self.allYList += y
        self.ax.plot(x, y, label=label, ls=ls, color=color)

    def finalAndSave(self, filename: str, xmax: float=None, ymin: float=None, ymax: float=None) -> None:
        # Core limits
        if (len(self.allYList) > 0):
            minValue = min(self.allYList) if ymin == None else ymin
            maxValue = max(self.allYList) if ymax == None else ymax
            avgValue = (maxValue+minValue)/2

            if (xmax == None):
                self.ax.vlines(
                    x=[rCore, roReflector, roReflectorSystem],
                    ymin=0 if self.isYlimZero else minValue,
                    ymax=maxValue,
                    ls='--', color='grey'
                )

                self.ax.text(
                    x=rCore/2, y=minValue*1.02,
                    s="Core", ha='center'
                )
                self.ax.text(
                    x=(roReflector+rCore)/2, y=minValue*1.02,
                    s="Relector", ha='center'
                )
                self.ax.text(
                    x=(roReflectorSystem+roReflector)/2, y=minValue*1.02,
                    s="Reflector\nsystem", ha='center'
                )

            self.ax.legend()

        if (self.isYlimZero):
            self.ax.set_ylim(0)
        elif (ymin != None and ymax != None):
            self.ax.set_ylim(ymin, ymax)
        else:
            self.ax.set_ylim(minValue, maxValue)


        if (xmax != None):
            self.ax.set_xlim(xmin=0, xmax=xmax)
        else:
            self.ax.set_xlim(xmin=0)


        self.fig.tight_layout()

        self.fig.savefig(filename)


#==============================================================================*
# Extract Results and Plot
#==============================================================================*

# Create figures and axis
figPower = Figure(
    "KIWI-B-4E core mid-plane radial relative fission distribution",
    r"Relative fission [$-$]",
)
figFlux = Figure(
    "KIWI-B-4E core mid-plane radial flux distribution",
    r"Neutron flux [n/cm$^2$/s]",
    isYlog=True
)
figTemperature = Figure(
    "KIWI-B-4E core mid-plane radial temperature distribution",
    r"Temperature [K]",
    isYlimZero=True
)


linestyles = ['-', '--', '-.', ':']

for filename, ls in zip(sys.argv[1:], linestyles):
    # File name
    folder = filename.split('/')[0]
    timeStep = filename.split('/')[-2]
    label = filename.split('/')[-1].split('_')[0]
    print(f"Process {folder} / {label} at time {timeStep}...")
    label = "GeN-Foam"

    # Update label if multiple files provided
    if (len(sys.argv[1:]) >= 2):
        label += " t="+timeStep

    # Extract all the parameter names from the filename (suppose OpenFOAM
    # autogeneration)
    names = ['r']+filename.split('/')[-1].split('.xy')[0].split('_')[1:]

    # Create empty list for each parameters
    data = {name: [] for name in names}

    # Open the file
    with open(filename, 'r') as file:
        # Extract all the data line by line
        for line in file.readlines():
            # Extract and format
            line = [float(e) for e in line.split()]

            # Fill in lists
            for i, name in enumerate(names):
                data[name].append(line[i])

        rRange = data['r']

        # Rescale variables
        for d in data:
            if ("flux" in d):
                # Convert from n/m2/s to n/cm2/s
                data[d] = listModifier(data[d], scale=1e-4)
                figFlux.plot(rRange, data[d], label=label+" "+d, ls=ls)

            elif ("power" in d):
                firstNon0idx = 0
                for i, d_ in enumerate(data[d]):
                    if (d_ > 0):
                        firstNon0idx = i
                        break
                data[d] = listModifier(data[d], scale=1/data[d][firstNon0idx])

                figPower.plot(rRange, data[d], label=label, ls=ls)

            elif ("T" == d[0]):
                figTemperature.plot(rRange, data[d], label=label+" "+d.split(".")[0], ls=ls)



# Plot reference values
plotRefFission(figPower.ax)

# Add vertical lines, legend and save
figPower.finalAndSave("fig_results_radialPowerDistribution.png", xmax=rCore, ymin=0.7, ymax=1.1)
figFlux.finalAndSave("fig_results_radialFluxDistribution.png")
figTemperature.finalAndSave("fig_results_radialTemperatureDistribution.png", xmax=rCore)


# plt.show()

#==============================================================================*
