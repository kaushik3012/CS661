### Import the required modules
###############################
from argparse import ArgumentParser
from vtk import *

### Parse the command line arguments
####################################
parser = ArgumentParser(description='Script to execute Isocontour Extraction on a 2D data')
parser.add_argument('--isovalue','-v',dest='isoval', help='Input the isovalue for contour extraction', required=True, type=float)
args = parser.parse_args()
isoval = args.isoval

### Create the reader for the data
##################################
reader = vtkXMLImageDataReader()            # create a reader
reader.SetFileName('Data/Isabel_2D.vti')    # set the file name
reader.Update()                             # update the reader
data = reader.GetOutput()                   # get the output of the reader

### Get the pressure array
##############################
dataArr = data.GetPointData().GetArray('Pressure')  # get the pressure array
numCells = data.GetNumberOfCells()                 # get the number of cells
cells = vtkCellArray()  # create a cell array
points = vtkPoints()    # create a point array

### Marching Squares Algorithm for contour extraction 
######################################################
for i in range(numCells):   # loop over all the cells

    cell = data.GetCell(i)  # get the cell with cellId

    ### Get the point ids of the cell
    #################################
    pid = [
            cell.GetPointId(0),
            cell.GetPointId(1),
            cell.GetPointId(3),
            cell.GetPointId(2)
        ]

    ### Get the scalar values of the points
    #######################################
    val = [
            dataArr.GetTuple1(pid[0]),
            dataArr.GetTuple1(pid[1]),
            dataArr.GetTuple1(pid[2]),
            dataArr.GetTuple1(pid[3])
        ]

    ### Get the 3D coordinates of the points
    ########################################
    p3d = [ 
            data.GetPoint(pid[0]),
            data.GetPoint(pid[1]),
            data.GetPoint(pid[2]),
            data.GetPoint(pid[3])
        ]        

    ### Check if the cell is completely inside or outside the isosurface
    ####################################################################
    if val[0]<isoval and val[1]<isoval and val[2]<isoval and val[3]<isoval:

        # if all the values are less than the isovalue, the cell is completely
        # outside the isosurface, so we can skip this cell
        continue

    elif val[0]>isoval and val[1]>isoval and val[2]>isoval and val[3]>isoval:

        # if all the values are greater than the isovalue, the cell is completely
        # inside the isosurface, so we can skip this cell
        continue

    else:

        # if the cell is partially inside the isosurface, we need to find the
        # intersection points and create a polyline for the cell 
        # (i.e. a line connecting the intersection points)

        c = 0         # counter for the number of intersection points

        ### Find the intersection points
        #################################
        if ((val[0]<=isoval and val[1]>isoval) or (val[0]>=isoval and val[1]<isoval)):
            px = ((val[0]-isoval)/(val[0]-val[1])) * (p3d[1][0]-p3d[0][0]) + p3d[0][0]
            py = ((val[0]-isoval)/(val[0]-val[1])) * (p3d[1][1]-p3d[0][1]) + p3d[0][1]
            p = (px,py,25)  
            points.InsertNextPoint(p)   # insert the point in the point array
            c+=1
        if ((val[1]<=isoval and val[2]>isoval) or (val[1]>=isoval and val[2]<isoval)):
            px = ((val[1]-isoval)/(val[1]-val[2])) * (p3d[2][0] - p3d[1][0]) + p3d[1][0]
            py = ((val[1]-isoval)/(val[1]-val[2])) * (p3d[2][1] - p3d[1][1]) + p3d[1][1]
            p = (px,py,25)
            points.InsertNextPoint(p)   # insert the point in the point array
            c+=1
        if ((val[2]<=isoval and val[3]>isoval) or (val[2]>=isoval and val[3]<isoval)):
            px = ((val[2]-isoval)/(val[2]-val[3])) * (p3d[3][0]-p3d[2][0]) + p3d[2][0]
            py = ((val[2]-isoval)/(val[2]-val[3])) * (p3d[3][1]-p3d[2][1]) + p3d[2][1]
            p = (px,py,25)
            points.InsertNextPoint(p)   # insert the point in the point array
            c+=1
        if ((val[3]<=isoval and val[0]>isoval) or (val[3]>=isoval and val[0]<isoval)):
            px = ((val[3]-isoval)/(val[3]-val[0])) * (p3d[0][0]-p3d[3][0]) + p3d[3][0]
            py = ((val[3]-isoval)/(val[3]-val[0])) * (p3d[0][1]-p3d[3][1]) + p3d[3][1]
            p = (px,py,25)
            points.InsertNextPoint(p)   # insert the point in the point array   
            c+=1    

        numPoints = points.GetNumberOfPoints() # get the number of points in the point array
        
        ### Create and Insert poly lines in the cell array
        ###################################################
        if c==2: 

            # If there are two intersection points, create a line
            # between the two pointS and insert the line in the cell array
            # Note: the point ids are the last two points in the point array
            #       since we have just inserted the two points in the previous step
            #       (numPoints-2) and (numPoints-1) are the ids of the last two 
            #       points in the point array

            polyLine = vtkPolyLine()    
            polyLine.GetPointIds().SetNumberOfIds(2)    
            polyLine.GetPointIds().SetId(0, numPoints-2)
            polyLine.GetPointIds().SetId(1, numPoints-1)
            cells.InsertNextCell(polyLine)

        elif c==4:

            # If there are four intersection points, create two lines
            # between the four points and insert the two lines in the cell array
            # Note: the point ids are the last four points in the point array
            #       since we have just inserted the four points in the previous step
            #       (numPoints-4), (numPoints-3), (numPoints-2) and (numPoints-1)
            #       are the ids of the last four points in the point array
            #       The two lines are created by connecting the first and second
            #       points and the third and fourth points

            polyLine1 = vtkPolyLine()
            polyLine1.GetPointIds().SetNumberOfIds(2)
            polyLine1.GetPointIds().SetId(0, numPoints-4)
            polyLine1.GetPointIds().SetId(1, numPoints-3)
            cells.InsertNextCell(polyLine1)

            polyLine2 = vtkPolyLine()
            polyLine2.GetPointIds().SetNumberOfIds(2)
            polyLine2.GetPointIds().SetId(0, numPoints-2)
            polyLine2.GetPointIds().SetId(1, numPoints-1)
            cells.InsertNextCell(polyLine2)

### Create a polydata
######################
pdata = vtkPolyData()

### Add points and cells to polydata
####################################
pdata.SetPoints(points)
pdata.SetLines(cells)

### Store the polydata into a vtkpolydata file with extension .vtp
###################################################################
writer = vtkXMLPolyDataWriter()
writer.SetInputData(pdata)
writer.SetFileName('xcontour.vtp')
writer.Write()