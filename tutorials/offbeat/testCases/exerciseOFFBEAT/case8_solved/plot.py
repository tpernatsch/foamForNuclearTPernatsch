import matplotlib.pyplot as plt

def read_data(file_path, data_columns, xAxis_factor=1, data_factor=1):
    """
    Reads xAxis and data values from a file, skipping lines starting with '#'.
    If data_columns is a single number, returns a single list of data.
    If data_columns is a list, returns a list of lists, each containing data for one column.
    Removes parentheses from vector data.
    """
    xAxis = []
    is_single_column = not isinstance(data_columns, list)
    if is_single_column:
        data_columns = [data_columns]
    data = [[] for _ in data_columns]  # Initialize a list of lists for data

    with open(file_path) as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split()
                xAxis.append(float(parts[0]) * xAxis_factor)

                for i, col in enumerate(data_columns):
                    # Remove parentheses and convert to float
                    value = parts[col].strip('()')
                    data[i].append(float(value) * data_factor)

    if is_single_column:
        return xAxis, data[0]  # Return a single list if only one column is provided
    else:
        return xAxis, data  # Return a list of lists if multiple columns are provided


def plot_data(ax, xAxis, data, axisLabel, legendLabel, color, line_style='-', 
              axis_id='y1', secondary_axes=None, y_limits=None, x_limits=None):
    """
    Plots data on the given axis, handling multiple secondary axes and axis limits.
    Includes the ability to specify the line style.
    Uses axis_id as a unique identifier for secondary axes.
    The first secondary axis ('y2') is attached to the main plot frame, 
    and subsequent secondary axes are pushed to the right.
    """
    if secondary_axes is None:
        secondary_axes = {}

    if axis_id == 'y1':
        ax.plot(xAxis, data, label=legendLabel, color=color, linestyle=line_style)
        ax.set_ylabel(axisLabel, color=color)
        ax.tick_params(axis='y', labelcolor=color)  # Set tick color for primary y-axis
        if y_limits:
            ax.set_ylim(y_limits)
        if x_limits:
            ax.set_xlim(x_limits)
    else:
        if axis_id not in secondary_axes:
            new_ax = ax.twinx()
            secondary_axes[axis_id] = new_ax
            new_ax.set_ylabel(axisLabel, color=color)

            # Adjust position for additional secondary axes beyond the first one
            if axis_id != 'y2':
                ax_position = ax.get_position()
                new_ax.set_position([ax_position.x0, ax_position.y0, 
                                     ax_position.width, ax_position.height])
                offset = 60 * (len(secondary_axes) - 1)  # Offset for additional axes
                new_ax.spines['right'].set_position(('outward', offset))

        existing_ax = secondary_axes[axis_id]
        existing_ax.plot(xAxis, data, label=legendLabel, color=color, linestyle=line_style)
        existing_ax.tick_params(axis='y', labelcolor=color)
        if y_limits:
            existing_ax.set_ylim(y_limits)
        if x_limits:
            existing_ax.set_xlim(x_limits)

    return ax, secondary_axes



# Paths
fuelOuter_path = "postProcessing/fuelOuter/0/"
fuelCentreline_path = "postProcessing/fuelCenterline/0/"
fuelTop_path = "postProcessing/fuelTop/0/"
cladTop_path = "postProcessing/cladTop/0/"

# Paths for radial profile data at specific times
radial_profile_paths = ["postProcessing/fuelRadialProfile/7.4/",
                        "postProcessing/fuelRadialProfile/168.4/",
                        "postProcessing/fuelRadialProfile/364/"]

# Positions
positions = ['bottom', 'mid', 'top']

# Constants for unit conversion
time_factor = 1  # Convert seconds to days
radialPos_factor = 1000 # Convert to mm

burnup_factor = 1.0 / 1000       # Convert to MWd/kg
gapWidth_factor = 1000000  # Convert to um
interfaceP_factor = 1.0 / 1000000  # Convert to MPa
hGap_factor = 1.0 / 1000000  # Convert to MW/m2/K
displacement_factor = 1000  # Convert to mm
stress_factor = 1.0 / 1e6  # Convert to mm

# Read fuel centerline temperature (FCT)
time, FCT = read_data(fuelCentreline_path + 'T', [1,2,3], time_factor)

# Read fuel centerline burnup
_, burnup = read_data(fuelCentreline_path + 'Bu', [1,2,3], time_factor, 
                        data_factor=burnup_factor)

# Read fuel outer gapWidth
_, gapWidth = read_data(fuelOuter_path + 'gapWidth', [1,2,3], time_factor, 
                        data_factor=gapWidth_factor)

# Read fuel outer hGap
_, hGap = read_data(fuelOuter_path + 'hGap', [1,2,3], time_factor, 
                        data_factor=hGap_factor)

# Read fuel outer interfaceP
_, interfaceP = read_data(fuelOuter_path + 'interfaceP', [1,2,3], time_factor, 
                        data_factor=interfaceP_factor)

# Read fuel and cladding axial expansion
_, fuelTopD = read_data(fuelTop_path + 'surfaceFieldValue.dat', [3], time_factor, 
                        data_factor=displacement_factor)
_, cladTopD = read_data(cladTop_path + 'surfaceFieldValue.dat', [3], time_factor, 
                          data_factor=displacement_factor)


# Initialize list to store radial profiles
radialT_profiles = []
radialSigmaR_profiles = []
radialSigmaH_profiles = []

# Read fuel radial profiles
for idx, pos in enumerate(positions):
    radialT_profiles.append([])
    radialSigmaR_profiles.append([])
    radialSigmaH_profiles.append([])
    for path in radial_profile_paths:
        radial_position, radialT = read_data(path + pos + '_T.xy', 1, radialPos_factor)
        radial_position, radialSigmaR = read_data(path + pos + '_sigma.xy', 1, 
            radialPos_factor, data_factor=stress_factor)
        radial_position, radialSigmaH = read_data(path + pos + '_sigma.xy', 4, 
            radialPos_factor, data_factor=stress_factor)
        
        radialT_profiles[idx].append(radialT) 
        radialSigmaR_profiles[idx].append(radialSigmaR) 
        radialSigmaH_profiles[idx].append(radialSigmaH) 


#### Plot 1 - FCT, Bu, gapWidth, hGap and interfaceP
fig1, ax1 = plt.subplots(figsize=(10, 5))  # Adjust figure size as needed
ax1.set_xlabel('Time, days')
ax1.set_title('T, Bu, Gap Width, Gap Conductance and Contact Pressure vs Time')
secondary_axes_1 = {}

# Plot FCTs on a primary axis
ax1, _ = plot_data(ax1, time, FCT[0], 'Temperature, K', 'Probe 1', 'tab:red', 
                   '-', 'y1', y_limits=[250, 2250])
ax1, _ = plot_data(ax1, time, FCT[1], 'Temperature, K', 'Probe 2', 'tab:red', 
                   '--', 'y1', y_limits=[250, 2250])
ax1, _ = plot_data(ax1, time, FCT[2], 'Temperature, K', 'Probe 3', 'tab:red', 
                   '-.', 'y1', y_limits=[250, 2250])

# Plot Bus  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, burnup[0], 'Burnup, MWd/kg', 
                   'Probe 1', 'tab:blue','-', 'y2', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, burnup[1], 'Burnup, MWd/kg', 
                   'Probe 2', 'tab:blue','--', 'y2', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, burnup[2], 'Burnup, MWd/kg', 
                   'Probe 3', 'tab:blue','-.', 'y2', secondary_axes_1)

# Plot gapWidth  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, gapWidth[0], 'Gap Width, $\mu$m', 
                   'Probe 1', 'tab:green','-', 'y3', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, gapWidth[1], 'Gap Width, $\mu$m', 
                   'Probe 2', 'tab:green','--', 'y3', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, gapWidth[2], 'Gap Width, $\mu$m', 
                   'Probe 3', 'tab:green','-.', 'y3', secondary_axes_1)

# Plot hGap  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, hGap[0], 'hGap, W/m$^2$/s', 
                   'Probe 1', 'tab:brown','-', 'y4', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, hGap[1], 'hGap, W/m$^2$/s', 
                   'Probe 2', 'tab:brown','--', 'y4', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, hGap[2], 'hGap, W/m$^2$/s', 
                   'Probe 3', 'tab:brown','-.', 'y4', secondary_axes_1)

# Plot interfaceP  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, interfaceP[0], 'Contact pressure, MPa', 
                   'Probe 1', 'tab:purple','-', 'y5', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, interfaceP[1], 'Contact pressure, MPa', 
                   'Probe 2', 'tab:purple','--', 'y5', secondary_axes_1)
ax1, secondary_axes_1 = plot_data(ax1, time, interfaceP[2], 'Contact pressure, MPa', 
                   'Probe 3', 'tab:purple','-.', 'y5', secondary_axes_1)

# Adjust subplots and layout
plt.subplots_adjust(right=0.75)  # Adjust the right boundary of the plot

# Adding legend for line types
line_styles = ['-', '--', '-.']
labels = ['Probe 1', 'Probe 2', 'Probe 3']
lines = [plt.Line2D([0], [0], color='black', linewidth=2, linestyle=ls)
         for ls in line_styles]
fig1.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.1, 0.9))
fig1.tight_layout()


#### Plot 2 - Radial Temperature Profiles (only at 364)
fig2, ax2 = plt.subplots()
ax2.set_xlabel('Radial Position, mm')
ax2.set_title('Radial Temperature Profile at 364 Days at Different Axial Levels')
secondary_axes_2 = {}

# Marker types for different times
markers = ['-', '--', '-.']  # Example markers: circle, square, triangle

# Loop over radial temperature profiles and plot them with different markers
for idx, pos in enumerate(positions):
    marker = markers[idx]
    time_label = '364'
    ax2, _ = plot_data(ax2, radial_position, radialT_profiles[idx][2], 'Temperature, K', 
        time_label, 'tab:red', marker, 'y1')


# Adding legend for line types
line_styles = ['-', '--', '-.']
labels = ['z = 1cm', 'z = 5cm', 'z = 10cm']
lines = [plt.Line2D([0], [0], color='tab:red', linewidth=2, linestyle=ls)
         for ls in line_styles]
fig2.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.15, 0.9))

fig2.tight_layout()


#### Plot 3 - Radial Stress Profiles (only at 364)
fig3, ax3 = plt.subplots()
ax3.set_xlabel('Radial Position, mm')
ax3.set_title('Radial Stress Profiles at 364 Days at Different Axial Levels')
secondary_axes_3 = {}

# Loop over radial stress radial profiles and plot them with different markers
for idx, pos in enumerate(positions):
    marker = markers[idx]
    time_label = '364'
    ax3, _ = plot_data(ax3, radial_position, radialSigmaR_profiles[idx][2], 'Radial Stress, MPa', 
                       time_label, 'k', marker, 'y1')   

# Adding legend for line types
line_styles = ['-', '--', '-.']
labels = ['z = 1cm', 'z = 5cm', 'z = 10cm']
lines = [plt.Line2D([0], [0], color='black', linewidth=2, linestyle=ls)
         for ls in line_styles]
fig3.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.15, 0.9))
fig3.tight_layout()


#### Plot 4 - Hoop Stress Profiles (only at 364)
fig3, ax3 = plt.subplots()
ax3.set_xlabel('Radial Position, mm')
ax3.set_title('Hoop Stress Profiles at 364 Days at Different Axial Levels')
secondary_axes_3 = {}

# Loop over radial stress radial profiles and plot them with different markers
for idx, pos in enumerate(positions):
    marker = markers[idx]
    time_label = '364'
    ax3, _ = plot_data(ax3, radial_position, radialSigmaH_profiles[idx][2], 'Hoop Stress, MPa', 
                       time_label, 'tab:cyan', marker, 'y1')   

# Adding legend for line types
line_styles = ['-', '--', '-.']
labels = ['z = 1cm', 'z = 5cm', 'z = 10cm']
lines = [plt.Line2D([0], [0], color='cyan', linewidth=2, linestyle=ls)
         for ls in line_styles]
fig3.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.15, 0.9))
fig3.tight_layout()


#### Plot 5 - Fuel and cladding axial elongation
fig5, ax5 = plt.subplots()
ax5.set_xlabel('Time, days')

# Plot gapWidth on a primary axis
ax5, _ = plot_data(ax5, time, fuelTopD[0], 'Displacement, mm', 'Fuel', 'tab:red', 
                   '-', 'y1')
ax5, _ = plot_data(ax5, time, cladTopD[0], 'Displacement, mm', 'Clad', 'k', 
                   '-', 'y1')

# Add the legend
ax5.legend(loc='upper left', bbox_to_anchor=(0.01, 0.99))
fig5.tight_layout()

plt.show()