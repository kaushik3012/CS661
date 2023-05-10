### Import VTK and other required modules
#########################################
from argparse import ArgumentParser
from vtk import *

### Parse the command line arguments
####################################
parser = ArgumentParser(description='Volume Rendering with Phong Shading')
parser.add_argument('--shading','-s',type=int, help='Toggle Phong shading', choices=[0,1], default=0)
args = parser.parse_args()

### Create the reader for the data
##################################
reader = vtkXMLImageDataReader()            # create a reader
reader.SetFileName('Data/Isabel_3D.vti')    # set the file name
reader.Update()                             # update the reader
data = reader.GetOutput()                   # get the output of the reader

### Create transfer mapping scalar value to opacity
###################################################
opacity_function = vtkPiecewiseFunction()
opacity_function.AddPoint(-4931.54,   1.0)
opacity_function.AddPoint(101.815, 0.002)
opacity_function.AddPoint(2594.97, 0.0)

### Create transfer mapping scalar value to color
#################################################
color_function = vtkColorTransferFunction()
color_function.AddRGBPoint(-4931.54, 0.0, 1.0, 1.0)
color_function.AddRGBPoint(-2508.95, 0.0, 0.0, 1.0)
color_function.AddRGBPoint(-1873.9, 0.0, 0.0, 0.5)
color_function.AddRGBPoint(-1027.16, 1.0, 0.0, 0.0)
color_function.AddRGBPoint(-298.031, 1.0, 0.4, 0.0)
color_function.AddRGBPoint(2594.97, 1.0, 1.0, 0.0)

### The property describes how the data will look
################################################
volume_property = vtkVolumeProperty()       
volume_property.SetColor(color_function)            # set the color function
volume_property.SetScalarOpacity(opacity_function)  # set the opacity function
volume_property.SetAmbient(0.5)                     # set the ambient property
volume_property.SetDiffuse(0.5)                     # set the diffuse property
volume_property.SetSpecular(0.5)                    # set the specular property
if args.shading:             # toggle phong shading
    volume_property.ShadeOn()   
else:   
    volume_property.ShadeOff()  
volume_property.SetInterpolationTypeToLinear()      # set the interpolation type

### The mapper / ray cast function know how to render the data
###############################################################
volume_mapper = vtkSmartVolumeMapper()
volume_mapper.SetInputData(data)    

### The actor object for the volume rendering
###############################################
volume_actor = vtkVolume()
volume_actor.SetMapper(volume_mapper)
volume_actor.SetProperty(volume_property)

### Create outline of the volume
################################
outlineFilter = vtkOutlineFilter()
outlineFilter.SetInputData(data)
outlineFilter.Update()

### Create mapper and actor for the outline
###########################################
outline_mapper = vtkPolyDataMapper() 
outline_mapper.SetInputConnection(outlineFilter.GetOutputPort())
outline_property = vtkProperty()
outline_property.SetColor(0,0,0)
outline_actor = vtkActor()
outline_actor.SetMapper(outline_mapper)   
outline_actor.SetProperty(outline_property) 

### Setup render window, renderer, and interactor
##################################################
renderer = vtkRenderer()
renderer.SetBackground(1,1,1)
renderWindow = vtkRenderWindow()
renderWindow.SetSize(1000,1000)
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderer.AddActor(volume_actor)
renderer.AddActor(outline_actor)
renderer.ResetCamera()

### Finally render the object
#############################
renderWindow.Render()
renderWindowInteractor.Start()