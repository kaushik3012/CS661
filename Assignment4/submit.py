import numpy as np
import vtk
from scipy.interpolate import griddata
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from argparse import ArgumentParser

def SRS(data, percent):
    # Convert the data into a numpy array
    data_array = np.array(data.GetPointData().GetScalars())
    dims = data.GetDimensions()
    data_array = data_array.reshape(dims[2], dims[1], dims[0])
    data_array = data_array.transpose(2,1,0)

    # Get the number of points in the dataset
    num_points = data.GetNumberOfPoints()

    # Calculate the number of points to sample
    num_sampled_points = int(num_points * percent / 100)

    # Randomly sample the points
    sampled_indices = np.random.choice(num_points, num_sampled_points, replace=False)
    eight_corner_points = [0, 
                            data.GetDimensions()[0]-1, 
                            data.GetNumberOfPoints()-1, 
                            data.GetNumberOfPoints()-data.GetDimensions()[0],
                            data.GetDimensions()[0]*(data.GetDimensions()[1]-1), 
                            data.GetDimensions()[0]*data.GetDimensions()[1]-1, 
                            data.GetNumberOfPoints()-data.GetDimensions()[0]*(data.GetDimensions()[1]-1)-1, 
                            data.GetNumberOfPoints()-data.GetDimensions()[0]*data.GetDimensions()[1]
                            ]
    
    # Add the eight corner points to the sampled indices
    sampled_indices = np.union1d(sampled_indices, eight_corner_points)

    # Sort the sampled indices in ascending order
    sampled_indices = np.sort(sampled_indices)

    # Create a new dataset for the sampled points
    sampled_data = vtk.vtkPolyData()
    points = vtk.vtkPoints()

    # Create an array to store the data values for each sampled point
    data_values = vtk.vtkFloatArray()
    data_values.SetNumberOfComponents(1)
    data_values.SetName("data")

    for i in sampled_indices:
        # Get the coordinates of the sampled point
        point = data.GetPoint(i)
        # Add the point to the new dataset
        points.InsertNextPoint(point)
        # Add the data value for the sampled point to the data values array
        data_values.InsertNextValue(data_array[int(point[0]), int(point[1]), int(point[2])])

    # Set the points and data values for the new dataset
    sampled_data.SetPoints(points)
    sampled_data.GetPointData().SetScalars(data_values)

    # Write the sampled data to a VTKPolyData file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName("sampled_data.vtp")
    writer.SetInputData(sampled_data)
    writer.Write()

    return sampled_data

def reconstruct_volume(sampled_data, method='nearest', dims=(250, 250, 50)):
    # Get the point coordinates and data values from the sampled data
    points = sampled_data.GetPoints()
    data_values = sampled_data.GetPointData().GetScalars()
    
    # Convert the point coordinates and data values to numpy arrays
    point_array = vtk_to_numpy(points.GetData())
    data_array = vtk_to_numpy(data_values)

    # Define the grid coordinates for the reconstruction
    grid_x, grid_y, grid_z = np.mgrid[:dims[0], :dims[1], :dims[2]]
    grid_points = (grid_x, grid_y, grid_z)

    # Perform the interpolation to reconstruct the volume data
    if method == 'nearest':
        # Use nearest neighbor interpolation
        reconstructed_data = griddata(point_array, data_array, grid_points, method='nearest')
        reconstructed_data = reconstructed_data.transpose(2,1,0)
    elif method == 'linear':
        # Use linear interpolation
        reconstructed_data = griddata(point_array, data_array, grid_points, method='linear')
        reconstructed_data = reconstructed_data.transpose(2,1,0)
        
        # Replace NaN values with nearest neighbor values
        nan_indices = np.isnan(reconstructed_data)
        if np.any(nan_indices):
            non_nan_indices = np.logical_not(nan_indices)
            non_nan_points = grid_points[non_nan_indices]
            non_nan_data = reconstructed_data[non_nan_indices]
            reconstructed_data[nan_indices] = griddata(non_nan_points, non_nan_data, grid_points[nan_indices], method='nearest')

    # Convert the reconstructed data back to a VTKImageData object
    vtk_data = numpy_to_vtk(reconstructed_data.flatten(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_image = vtk.vtkImageData()
    vtk_image.SetDimensions(dims[0], dims[1], dims[2])
    vtk_image.GetPointData().SetScalars(vtk_data)

    # Write the reconstructed data to a VTI file
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(f"reconstructed_data_{method}.vti")
    writer.SetInputData(vtk_image)
    writer.Write()

    return vtk_image

def compute_SNR(arrgt, arr_recon):
        diff = arrgt - arr_recon
        sqd_max_diff = (np.max(arrgt) - np.min(arrgt))**2
        snr = 10*np.log10(sqd_max_diff/np.mean(diff**2))
        return snr

if __name__=="__main__":

    # Parse the command line arguments
    parser = ArgumentParser(description='Script to execute Sampling and Volume Reconstruction')
    parser.add_argument('--samp_percent','-s',dest='samp_percent', help='Input the percent of points required for Simple Random Sampling', required=True, type=float)
    parser.add_argument('--recon_method','-r',dest='recon_method', help='Input the reconstruction method', required=True, type=str)
    args = parser.parse_args()
    samp_percent = args.samp_percent
    recon_method = args.recon_method

    # Load the input data from a VTKPolyData file
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName("Data/Isabel_3D.vti")
    reader.Update()
    input_data = reader.GetOutput()

    # Sample the input data using simple random sampling with a sampling percentage of 10%
    output_data = SRS(input_data, samp_percent)
    print('Sampling complete')

    # Load the input data from a VTKPolyData file
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName("sampled_data.vtp")
    reader.Update()
    output_data = reader.GetOutput()

    reconstr_vol = reconstruct_volume(output_data, method=recon_method)
    print('Reconstruction complete')

    gt_np = vtk_to_numpy(input_data.GetPointData().GetScalars())
    recon_np = vtk_to_numpy(reconstr_vol.GetPointData().GetScalars())
    print(f'SNR: {compute_SNR(gt_np, recon_np)}')