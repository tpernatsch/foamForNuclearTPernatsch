import os
import subprocess
from PIL import Image

# Path to the "Frames" folder containing the PNG files
frames_folder = "Frames"

# Output directory for the MP4 file
output_directory = "Output"

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Get the list of PNG files in the "Frames" folder
png_files = sorted([f for f in os.listdir(frames_folder) if f.endswith(".png")])

# Calculate dimensions based on the aspect ratio of the first image
first_img_path = os.path.join(frames_folder, png_files[0])
first_img = Image.open(first_img_path)
aspect_ratio = first_img.width / first_img.height

# Define the target width and height for resizing
target_width = 3000  # Adjust this to your desired width
# target_height = int(target_width / aspect_ratio)
target_height = 1694

# Convert PNG frames to MP4 using ffmpeg with proper resizing
mp4_filepath = os.path.join(output_directory, "ramp.mp4")
input_pattern = os.path.join(frames_folder, "ramp.%04d.png")

# The following command uses ffmpeg to convert the PNG frames to MP4 format.
# Adjust the `-r` option to control the frame rate (frames per second).
ffmpeg_command = (
    f"ffmpeg -framerate 20 -i {input_pattern} "
    f"-vf scale={target_width}:{target_height} "
    f"-c:v libx264 -pix_fmt yuv420p {mp4_filepath}"
)
subprocess.run(ffmpeg_command, shell=True)

print(f"MP4 animation created and saved as '{mp4_filepath}'.")
