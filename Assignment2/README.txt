******************************************************************************************
    Assignment 2 - CS661A
    Isocontour and Volume Visualization
    Submitted by: Kaushik Raj V Nadar (200499)
******************************************************************************************

******************************************************************************************
   Question 1: 2D Isocontour Extraction
******************************************************************************************
The question 1 of the assignment is implemented as a python script in q1.py file. 
The script takes the isovalue for contouring as a command line argument. 
The script generates a .vtp file which can be opened in paraview to 
visualize the contour. By default, the name of the output file is 'xcontour.vtp'.

If the required contour has a isovalue of 100, then the script is run as follows:
python q1.py --isovalue=100.0
or
python q1.py -v 100.0
*******************************************************************************************

*******************************************************************************************
   Question 2: VTK Volume Rendering and Transfer Function
*******************************************************************************************
The question 2 of the assignment is implemented as a python script in q2.py file.
The script takes an input parameter (0 or 1) as a command line argument to toggle 
the Phong shading on and off.
The script renders a volume using the transfer function provided in the assignment.

If the Phong shading is to be turned off, then the script is run as follows:
python q2.py --shading=0
or
python q2.py -s 0

If the Phong shading is to be turned on, then the script is run as follows:
python q2.py --shading=1
or
python q2.py -s 1
*********************************************************************************************
