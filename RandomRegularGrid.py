"""Provides a scripting component.
    Inputs:
        P: The polygon to subdivide
        GS: Integer of the minimum grid size
        M: Maximun Cells in a Plot
        S: Random Seed
    Output:
        D: The Subdivisions"""

__author__ = "crubiogarcia"
__version__ = "2023.06.14"

import rhinoscriptsyntax as rs
import Rhino.Geometry as rg
import scriptcontext as sc
import random
import operator
import math
import time

#Create timer
class Timer(object):
    
    """ A simple profiler """
    
    def start(self):
        # Start the timer
        self.startTime = time.time()
        
    def stop(self):
        # Get elapsed time
        elapsedTime = time.time() - self.startTime
        elapsedMilliseconds = elapsedTime*1000
        print str(round(elapsedMilliseconds,3)) + "ms"
        return elapsedMilliseconds

#Main Function
def Main():
    
    global M
    M = int(M)
    global S
    S = int(S)
    global P
    P = P
    global GS
    GS = int(GS)
    
    #Chek if Polyline is closed
    if rs.IsCurveClosed(P) == False:
        print("Polygon must be closed")
        return 1
    
    #Chek if P is a Polyline
    if rs.IsPolyline(P) == False:
        print("P must be a polyline")
        return 2
    
    #Chek if P is Planar
    if rs.IsCurvePlanar(P) == False:
        print("P must be planar")
        return 3

    if GS <= 0:
        print("Grid Size must be larger than 0")
        return 4
        
    if M <= 0:
        print("M must be larger than 0")
        return 5
        
        

    #get curve lengths
    def curve_lengths(curves):
        return [rs.CurveLength(curve) for curve in curves]
        
    #draw rectangle
    def rectangle(corner_0, corner_1):
        return rs.AddRectangle(rs.WorldXYPlane(), corner_0, corner_1)
        
       
    #create class for the grid cells
    class cell:
        def __init__(self,brep):
            self.brep = brep
            self.center = rs.SurfaceAreaCentroid(brep)[0]
            self.edges = rs.DuplicateSurfaceBorder(brep)
            self.vertices = rs.PolylineVertices(self.edges)[:4]
            
    #Subdivide between two lines 
    def get_subdiv(lines, subdiv, lengths):
        #creates curves in between two lines witha desired minimum spacing
        num = min([lengths[0]/subdiv, lengths[0]/subdiv])
        a1 = rs.DivideCurve(lines[0],num,True,True)
        a2 = rs.DivideCurve(lines[1],num,True,True)
        a2.reverse()
        ln = [rs.AddLine(a1[f], a2[f]) for f in range(len(a1))]
        return ln
        
    #Distance between Point and List of Points
    def distances(point, point_list):
        dist = [rs.Distance(point, i) for i in point_list]
        return dist
    
    #GET GRID, NUMBER OF ROWS, COLUMNS
    def get_grid(closed_crv, cell_size):
        
        #BOUNDING RECTANGLE
        #get reference plane based on largest edge of polygon
        edges = rs.ExplodeCurves(closed_crv)
        length = curve_lengths(edges)
        
        idx = length.index(max(length))
        aux_plane = rs.PlaneFromNormal([0,0,0],[0,1,0])

        start = rs.CurveStartPoint(edges[idx])
        end = rs.CurveEndPoint(edges[idx])
        vect = end-start

        p_start = rs.PlaneClosestPoint(aux_plane, start)
        p_end = rs.PlaneClosestPoint(aux_plane, end)
        p_vect = p_end-p_start

        ang = math.degrees(math.atan2(p_vect[1], p_vect[0]) - math.atan2(vect[1], vect[0]))
        
        sr = rs.AddPlanarSrf(edges)

        obj = []
        cent = rs.SurfaceAreaCentroid(sr)[0]
        curve = rs.RotateObjects(closed_crv, cent, ang)
        edges = rs.RotateObjects(edges, cent, ang)
    
        vect = rs.VectorCreate(rs.CurveEndPoint(edges[idx]),rs.CurveStartPoint(edges[idx])) 
        aux_crv = rs.OffsetCurve(edges[idx], vect, 10)
        
        crv_pl = []
        crv_pl.append(edges[idx])
        crv_pl.append(aux_crv)
    
        pts = rs.SurfacePoints(rs.AddLoftSrf(crv_pl))
        #get the bounding box
        bb = rs.BoundingBox(closed_crv,rs.PlaneFromPoints(pts[0],pts[1],pts[2]))
        
        #get the bounding rectangle
        box = rs.AddPolyline(bb[:5])
        box_srf = rs.AddPlanarSrf(box)
    
        #GRID
        #get number of subdivisions
        lns = rs.ExplodeCurves(box)
        lines_a = [lns[0], lns[2]]
        lines_b = [lns[1], lns[3]]
        lng_a = curve_lengths(lines_a)
        lng_b = curve_lengths(lines_b)
    
        grid_a = get_subdiv(lines_a, cell_size, lng_a)
        grid_b = get_subdiv(lines_b, cell_size, lng_b)
        
        gr_size = [len(grid_a)-1, len(grid_b)-1]
        gr_size_mx_idx = gr_size.index(max(gr_size))
        gr_size_min_idx = gr_size.index(min(gr_size))
    
        grid = []
        for l in grid_a:
            grid.append(l)
        for l in grid_b:
            grid.append(l)
        
        aux_ln = rs.AddLine([0,0,0],[0,0,5])
    
        cutters=[rs.coercecurve(ID) for ID in grid]
        bre =rs.coercebrep(box_srf)
        to_cut = bre.Faces[0]
        
        sp_result = rg.BrepFace.Split(to_cut, cutters,sc.doc.ModelAbsoluteTolerance)
        split_result = sp_result.Faces
        
        cel = []
        for f in split_result:
            cel.append(f.DuplicateFace(False))
        
        cls = []
        for bp in cel: 
            cls.append(sc.doc.Objects.AddBrep(bp))
            
        grid = cls
    
        cs = []
        for c in range(len(cls)):
            ce = cell(cls[c])
            cs.append(ce)
    
        #sort cells in the grid
        cs = sorted(cs, key=lambda k: k.center[gr_size_mx_idx])
        
        nest = [cs[i:i+gr_size[gr_size_min_idx]] for i in range(0, len(cs), gr_size[gr_size_min_idx])]
        for l in nest:
            l = sorted(l, key=lambda k: k.center[gr_size_mx_idx])
        
        fl = []
        for t in range(len(nest)):
            nest[t] = sorted(nest[t], key=lambda k: k.center[gr_size_min_idx])
            for p in range(len(nest[t])):
                fl.append(nest[t][p])
    
        cells = []
        for c in fl:
            cells.append(c)
            
        return cells, gr_size, box, cent, ang, sr
        
    #GET NEIGHBOURS
    def neighbours(cell_coordinates):
        neighbours = []
        min_x, min_y = M - 1, M - 1
        max_x, max_y = GR_SIZE[0] - (1 + M), GR_SIZE[1] - (1 + M)
    
        for i in range(-M + 1, M):
            for j in range(-M + 1, M):
                x, y = cell_coordinates[0] + i, cell_coordinates[1] + j
                if 0 <= x <= GR_SIZE[0]-1 and 0 <= y <= GR_SIZE[1]-1:
                    neighbours.append([x, y])
    
        return neighbours
    
    def round_neighbours(cell_coordinates, neighbours_coordinates, step):
        n = []
        for i in neighbours_coordinates:
            if i[0] == cell_coordinates[0] + step or i[0] == cell_coordinates[0] - step:
                if cell_coordinates[1] - step <= i[1] and i[1] <= cell_coordinates[1] + step:
                    n.append(i)
            elif i[1] ==  cell_coordinates[1] + step or i[1] ==  cell_coordinates[1] - step:
                if cell_coordinates[0] - step <= i[0] and i[0] <= cell_coordinates[0] + step:
                    n.append(i)
        return n
    
    def check_conditions(list_1, list_2,operator_a, operator_b, equal_1, equal_2):
        if equal_2 == False:
            for i in range(len(list_2)):
                for j in list_1:
                    if operator_a(j[0], list_2[i][0]):
                        if equal_1 != True and operator_b(j[1], list_2[i][1]):
                            to_remove.append(j)
                        elif equal_1 == True:
                            to_remove.append(j)
        elif equal_2 == True:
            for i in range(len(list_2)):
                for j in list_1:
                    if operator_a(j[1], list_2[i][1]):
                        to_remove.append(j)
            
            
    
    
    #Set Grid 
    aux = get_grid(P, GS)
    cell_list = aux[0]
    #Make the Grid not Nested
    grid =[]
    for c in cell_list:
        grid.append(c.brep)
    
    #Set Grid Size
    GR_SIZE = aux[1]
    
    
    #Nest List sorted
    nested = [cell_list[i:i+min(GR_SIZE)] for i in range(0, len(cell_list), min(GR_SIZE))]
    qu = [j for i in nested for j in i]
    #Create a list of visited
    visited = [[False for j in range(len(i))] for i in nested]
    
    #Create a flat index list
    auxi = range(GR_SIZE[0]*GR_SIZE[1])
    #Pop random element of the index list
    
    rectangles = []
    s = S
    
    while len(auxi)>0:
        random.seed(s)
        ix = auxi.pop(auxi.index(random.choice(auxi)))
        #Translate index into the nested list indexes 
        a, b = ix//GR_SIZE[1], ix%GR_SIZE[1]
        #Get Cell Coordinates
        cord = [a,b]
        
        #Get Cell Neighbours
        nghbs = neighbours(cord)
        
        #CHECK IF ANY OF THE NEIGHBOURS ARE VISITED
        #Create list
        check = [False for i in nghbs]
        larger = {'larger': [], 'smaller': [], 'equal': []}
        smaller = {'larger': [], 'smaller': [], 'equal': []}
        equal = {'larger': [], 'smaller': []}
        
        #Create list of elements to be removed
        to_remove = []
        for i in range(M-1):
            #Get Closest Neighbours
            r_nghbs = round_neighbours(cord, nghbs, i+1)
            #Append to List if visited 
            for i in range(len(r_nghbs)):
                if visited[r_nghbs[i][0]][r_nghbs[i][1]] == True:
                    check[i] = True
                    if r_nghbs[i][0]>cord[0] and r_nghbs[i][1]>cord[1]:
                        larger['larger'].append(r_nghbs[i])
                    elif r_nghbs[i][0]>cord[0] and r_nghbs[i][1]<cord[1]:
                        larger['smaller'].append(r_nghbs[i])
                    elif r_nghbs[i][0]>cord[0] and r_nghbs[i][1]==cord[1]:
                        larger['equal'].append(r_nghbs[i])
                
                    elif r_nghbs[i][0]<cord[0] and r_nghbs[i][1]>cord[1]:
                        smaller['larger'].append(r_nghbs[i])
                    elif r_nghbs[i][0]<cord[0] and r_nghbs[i][1]<cord[1]:
                        smaller['smaller'].append(r_nghbs[i])
                    elif r_nghbs[i][0]<cord[0] and r_nghbs[i][1]==cord[1]:
                        smaller['equal'].append(r_nghbs[i])      
                
                    elif r_nghbs[i][0] == cord[0] and r_nghbs[i][1]>cord[1]:
                        equal['larger'].append(r_nghbs[i])
                    elif r_nghbs[i][0] == cord[0] and r_nghbs[i][1]<cord[1]:
                        equal['smaller'].append(r_nghbs[i])
            
            #Remove elements that dont meet the conditions as a potential corner
            for i in range(len(check)):
                if check[i] == True:
                    to_remove.append(i)
        
        to_remove = list(set(to_remove))
        
        check_conditions(nghbs, larger['larger'], operator.ge, operator.ge, False, False)
        check_conditions(nghbs, larger['smaller'], operator.ge, operator.le, False, False)
        check_conditions(nghbs, larger['equal'], operator.ge, operator.ge, True, False)
        
        check_conditions(nghbs, smaller['larger'], operator.le, operator.ge, False, False)
        check_conditions(nghbs, smaller['smaller'], operator.le, operator.le, False, False)
        check_conditions(nghbs, smaller['equal'], operator.le, operator.le, True, False)
        
        check_conditions(nghbs, equal['larger'], operator.ge, operator.ge, True, True)
        check_conditions(nghbs, equal['smaller'], operator.le, operator.le, True, True)
        
        not_visited = [i for i in nghbs if i not in to_remove]
        not_visited.append([a,b])
        
        random.seed(s+1)
        cc = random.choice(not_visited)
        n_c = nested[cc[0]][cc[1]]
        C = n_c.brep
        
        #draw rectangle
        if n_c == nested[cord[0]][cord[1]] :
            rect = rectangle(n_c.vertices[0], n_c.vertices[2])
            
        else:
            
            dis_0 = distances(nested[cord[0]][cord[1]].center, n_c.vertices)
            i_0 = dis_0.index(max(dis_0))
            corner_0 = n_c.vertices[i_0]
            
            dis_1 = distances(corner_0, nested[cord[0]][cord[1]].vertices)
            i_1 = dis_1.index(max(dis_1))
            corner_1= nested[cord[0]][cord[1]].vertices[i_1]
            
            rect = rectangle(corner_0, corner_1)
        
        rectangles.append(rect)
        
        #Mark as visited all cells inside rectangle
        vis_ind = []
        for i in nghbs:
            if rs.PointInPlanarClosedCurve(nested[i[0]][i[1]].center, rect) == True:
                visited[i[0]][i[1]] = True
                vis_ind.append(i[0]*GR_SIZE[1] + i[1])
        
        #Remove visited indexes from index list
        auxi = [i for i in auxi if i not in vis_ind]
        s = 3*s/2
    
    
    #Rotate geometry into initial position
    rs.RotateObjects(rectangles, aux[3], -aux[4])
    
    #SPLIT POLYGON WITH THE RECTANGLES
    rects = rs.AddPlanarSrf(rectangles)
    
    plots = []
    for i in rects:
        srfs = rs.SplitBrep(aux[5], i)
        if srfs:
            centroids = [rs.SurfaceAreaCentroid(j)[0] for j in srfs]
            for k in range(len(centroids)):
                if rs.PointInPlanarClosedCurve(centroids[k], P):
                    plots.append(rs.DuplicateSurfaceBorder(srfs[k]))
                    
    subdiv = [j for i in plots for j in i]
    
    return subdiv

#Set timer
t = Timer()
#Start timer
t.start()

#Run script
D = Main()

#Stop timer
t.stop()
