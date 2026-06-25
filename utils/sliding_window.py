import numpy as np
import mrcfile
import csv


def generate_subvolumes(
    input_file, output_prefix, subvolume_size=(20, 20, 20), step_size=(15, 15, 15)
):
    """
    Generate subvolumes from a 3D tomogram, cropping it into smaller 3D volumes using a sliding window approach.

    Parameters:
    ----------
    input_file : str
        Path to the input .mrc file containing the 3D tomogram.
    output_prefix : str
        Prefix for names of the output subvolume .mrc files.
    subvolume_size : tuple of int, optional
        Size of each subvolume in (z_dim, y_dim, x_dim). Default is (100, 100, 100).
    step_size : tuple of int, optional
        Step size (stride) for the sliding window in (z_dim, y_dim, x_dim). Default is (30, 30, 30).

    Returns:
    --------
    list of str
        List of file paths for all generated .mrc subvolumes.

    Side Effects:
    -------------
    - Saves extracted subvolumes to .mrc files named with the specified prefix and their sliding window's starting position.

    Example Usage:
    --------------
    # Example Input File: "your_volume_file.mrc"
    # Example Output Prefix: "subvolume"
    # Example subvolume size: (100, 100, 100)
    # Example step size: (30, 30, 30)

    files = generate_subvolumes("your_volume_file.mrc", "subvolume", (100, 100, 100), (30, 30, 30))
    """
    # Open the MRC file
    with mrcfile.open(
        input_file, mode="r+"
    ) as mrc:  # changing from above line, to mode='r+', read and write mode
        tomogram = mrc.data  # Load the 3D data
        # here save the coordinate w.r.t the parent tomogram - so as to save global location
        retained_coords = mrc.header.origin  # HERE - include origin retainer
        retained_coords = retained_coords.ravel()
        parent_x = retained_coords[0][0]  # first without int()
        parent_y = retained_coords[0][1]
        parent_z = retained_coords[0][2]
        ## could try ..
        # p_new = (0., 0., 0.)
        # mrc.header.origin = p_new# this resets original tomogram which is good but tests will still fail

        # Ensure the tomogram is 3D
        if tomogram.ndim != 3:
            raise ValueError("Input file does not contain a 3D tomogram.")

        # Get the dimensions of the original volume
        z_dim, y_dim, x_dim = tomogram.shape
        sub_z, sub_y, sub_x = subvolume_size
        step_z, step_y, step_x = step_size

        # --
        generated_files = []  # List to store output file paths
        csv_log = f"{subvolume_size}_{step_size}.csv"  # include also a csv output
        csv_header = [
            "file name",
            "x-start",
            "x-end",
            "y-start",
            "y-end",
            "z-start",
            "z-end",
        ]
        with open(csv_log, "w", newline="") as file:  # open here to establish header
            writer = csv.writer(file)
            writer.writerow(csv_header)
        # --

        # Sliding window over the tomogram
        for z_start in range(0, z_dim - sub_z + 1, step_z):
            for y_start in range(0, y_dim - sub_y + 1, step_y):
                for x_start in range(0, x_dim - sub_x + 1, step_x):
                    # Define the subvolume boundaries
                    z_end = z_start + sub_z
                    y_end = y_start + sub_y
                    x_end = x_start + sub_x

                    # Extract the subvolume
                    subvolume = tomogram[z_start:z_end, y_start:y_end, x_start:x_end]

                    # Define a unique name for the output file
                    output_file = f"{output_prefix}_x{x_start}-{x_end}_y{y_start}-{y_end}_z{z_start}-{z_end}.mrc"  # alter file name convention to include coordinate min/max
                    pertaining_row = []
                    pertaining_row.append(output_file)
                    pertaining_row.append(x_start)
                    pertaining_row.append(x_end)
                    pertaining_row.append(y_start)
                    pertaining_row.append(y_end)
                    pertaining_row.append(z_start)
                    pertaining_row.append(z_end)
                    with open(
                        csv_log, mode="a", newline=""
                    ) as file:  # open here to input file and coords
                        writer = csv.writer(file)
                        writer.writerow(pertaining_row)

                    # Save the subvolume to .mrc format
                    with mrcfile.new(output_file, overwrite=True) as mrc_out:
                        mrc_out.set_data(subvolume.astype(np.float32))
                        mrc_out.voxel_size = (
                            mrc.voxel_size
                        )  # Copy voxel size from input file for consistency
                        # mrc_out.header.origin = retained_coords# ensure overall position is saved
                        voxels = mrc.voxel_size
                        voxels = voxels.ravel()
                        voxel_val = voxels[0][0]
                        mrc_out.header.origin = (
                            parent_x + (x_start * voxel_val),
                            parent_y + (y_start * voxel_val),
                            parent_z + (z_start * voxel_val),
                        )
                        # add coordinate storage ... csv?

                    print(f"Generated subvolume: {output_file}")
                    generated_files.append(output_file)

        return generated_files


# Example usage:
# generate_subvolumes(
#     "Your_Tomogram_in_MRC_format_goes_here",
#     "sub_volume",
#     subvolume_size=(219, 219, 219),
#     step_size=(110, 110, 110),
# )
