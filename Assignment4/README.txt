******************************************************************************************
    Assignment 4 - CS661A
    Sampling and Reconstruction
    Submitted by: Kaushik Raj V Nadar (200499)
******************************************************************************************

The assignment is implemented as a python script in submit.py file. 
The script takes the sampling percentage and volume reconstruction method as command line arguments. 
The script generates a .vtp file for sampled data and .vti file for reconstructed data which can be opened in paraview to 
visualize.

If the required sampling percentage is 1% and reconstruction method is nearest, then the script is run as follows:
python submit.py --samp_percent=1 --recon_method=nearest
or
python submit.py -s=1 -r=nearest

If the required sampling percentage is 5% and reconstruction method is linear, then the script is run as follows:
python submit.py --samp_percent=5 --recon_method=linear
or
python submit.py -s=5 -r=linear

*******************************************************************************************

