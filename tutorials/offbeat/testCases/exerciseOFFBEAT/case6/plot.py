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
radial_profile_paths = ["postProcessing/fuelRadialProfile/7.04/",
                        "postProcessing/fuelRadialProfile/168.04/",
                        "postProcessing/fuelRadialProfile/364/"]

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
time, FCT = read_data(fuelCentreline_path + 'T', [1], time_factor)

# Read fuel centerline burnup
_, burnup = read_data(fuelCentreline_path + 'Bu', [1], time_factor, 
                        data_factor=burnup_factor)

# Read fuel outer gapWidth
_, gapWidth = read_data(fuelOuter_path + 'gapWidth', [1], time_factor, 
                        data_factor=gapWidth_factor)

# Read fuel outer hGap
_, hGap = read_data(fuelOuter_path + 'hGap', [1], time_factor, 
                        data_factor=hGap_factor)

# Read fuel outer interfaceP
_, interfaceP = read_data(fuelOuter_path + 'interfaceP', [1], time_factor, 
                        data_factor=interfaceP_factor)

# Initialize list to store radial profiles
radialT_profiles = []
radialSigmaR_profiles = []
radialSigmaH_profiles = []

# Read fuel radial profile
for path in radial_profile_paths:
    radial_position, radialT = read_data(path + 'line_T.xy', 1, radialPos_factor)
    radial_position, radialSigmaR = read_data(path + 'line_sigma.xy', 1, 
        radialPos_factor, data_factor=stress_factor)
    radial_position, radialSigmaH = read_data(path + 'line_sigma.xy', 4, 
        radialPos_factor, data_factor=stress_factor)
    
    radialT_profiles.append(radialT) 
    radialSigmaR_profiles.append(radialSigmaR) 
    radialSigmaH_profiles.append(radialSigmaH) 

#### Plot 1 - FCT, Bu, gapWidth, hGap and interfaceP
fig1, ax1 = plt.subplots(figsize=(10, 5))  # Adjust figure size as needed
ax1.set_xlabel('Time, days')
secondary_axes_1 = {}

# Plot FCTs on a primary axis
ax1, _ = plot_data(ax1, time, FCT[0], 'Temperature, K', 'FCT', 'tab:red', 
                   '-', 'y1', y_limits=[250, 2500])
# Plot Bu  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, burnup[0], 'Burnup, MWd/kg', 
                   'Burnup', 'tab:blue','-', 'y2', secondary_axes_1)
# Plot gapWidth  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, gapWidth[0], 'Gap Width, $\mu$m', 
                   'Gap Width', 'tab:green','-', 'y3', secondary_axes_1)
# Plot gapWidth  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, hGap[0], 'hGap, W/m$^2$/s', 
                   'hGap', 'tab:brown','-', 'y4', secondary_axes_1)
# Plot gapWidth  on a secondary axis
ax1, secondary_axes_1 = plot_data(ax1, time, interfaceP[0], 'Contact pressure, MPa', 
                   'Contact pressure', 'tab:purple','-', 'y5', secondary_axes_1)

# Adjust subplots and layout
plt.subplots_adjust(right=0.75)  # Adjust the right boundary of the plot
fig1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig1.tight_layout()

#### Plot 2 - Radial Temperature Profiles
fig2, ax2 = plt.subplots()
ax2.set_xlabel('Radial Position, mm')
ax2.set_title('Radial Temperature Profile at Different Times')
secondary_axes_2 = {}

# Marker types for different times
markers = ['-', '--', '-.']  # Example markers: circle, square, triangle

# Loop over radial temperature profiles and plot them with different markers
for radialT, marker, path in zip(radialT_profiles, markers, radial_profile_paths):
    time_label = path.split('/')[-2]  # Extract time label from the path
    ax2, _ = plot_data(ax2, radial_position, radialT, 'Temperature, K', 
        time_label, 'tab:red', marker, 'y1')

fig2.legend(loc='upper left', bbox_to_anchor=(0.15, 0.95))
fig2.tight_layout()

#### Plot 3 - Radial Stress Profiles
fig3, ax3 = plt.subplots()
ax3.set_xlabel('Radial Position, mm')
ax3.set_title('Radial Stress Profiles at Different Times')
secondary_axes_3 = {}

# Loop over radial stress radial profiles and plot them with different markers
for radialSigmaR, marker, path in zip(radialSigmaR_profiles, markers, radial_profile_paths):
    time_label = path.split('/')[-2]  # Extract time label from the path
    ax3, _ = plot_data(ax3, radial_position, radialSigmaR, 'Radial Stress, MPa', 
                       time_label, 'k', marker, 'y1')

fig3.legend(loc='upper left', bbox_to_anchor=(0.15, 0.95))
fig3.tight_layout()

#### Plot 4 - Hoop Stress Profiles
fig4, ax4 = plt.subplots()
ax4.set_xlabel('Radial Position, mm')
ax4.set_title('Hoop Stress Profiles at Different Times')
secondary_axes_4 = {}

# Loop over hoop stress radial profiles and plot them with different markers
for radialSigmaH, marker, path in zip(radialSigmaH_profiles, markers, radial_profile_paths):
    time_label = path.split('/')[-2]  # Extract time label from the path
    ax4, _ = plot_data(ax4, radial_position, radialSigmaH, 
        'Hoop Stress, MPa', time_label, 'tab:cyan', marker, 'y1')

fig4.legend(loc='upper left', bbox_to_anchor=(0.15, 0.95))
fig4.tight_layout()

plt.show()