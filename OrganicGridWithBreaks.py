"""Provides a scripting component.
    Inputs:
        P: Closed polyline
        M: Minimum length of polygons
        I: Number of iterations
        S: Random seed
    Output:
        a: The a output variable"""

__author__ = "crubiogarcia"
__version__ = "2023.07.03"

import rhinoscriptsyntax as rs
import random
import ghpythonlib.treehelpers as gt

def check_curve_orientation(tt):
    if rs.ClosedCurveOrientation(tt) != -1:
        rs.ReverseCurve(P)
        return

#Function to subdivide closed curve
def subd(crvs, mn, sd):
    
    #Get breps from closed planar surfaces
    breps = rs.AddPlanarSrf(crvs)
    
    #Get brep edges
    bs = [rs.DuplicateEdgeCurves(brep) for brep in breps]

    sides = [ls for ls in bs]

    #Get edges lengths
    lengths = []
    for ls in sides:
        aux = [rs.CurveLength(pl) for pl in ls]
        lengths.append(aux)
    
    #Get max length edge
    mx = [max(m) for m in lengths]
    
    idx = [lengths[m].index(mx[m]) for m in range(len(lengths))]
            
    nn = []

    for j in range(len(mx)):
        if mx[j]>mn*2 and mx[j] != mn:
            #Get normalized min and max parameters in curve
            st = mn/mx[j]
            ed = 1-mn/mx[j]

            #Get a random choice parameter inside the range established
            random.seed(sd)
            repa = round(random.uniform(st, ed), 2)
            #Get from normalized parameter the parameter not normalized
            param = rs.CurveParameter(sides[j][idx[j]], repa)

            #Get point in curve
            pt = rs.EvaluateCurve(sides[j][idx[j]], param)

            #Get normal vector in curve
            tangent = rs.CurveTangent(sides[j][idx[j]], param)
            norm = rs.VectorRotate(tangent, 90, [0,0,1])
            
            #Draw line to split elements with normal direction
            line = rs.AddLine(pt, pt+norm)
            long = rs.ExtendCurveLength(line, 0, 2, max(mx)*10000)
            
            #getting points in the curve to rebuild the curve just inside the curve
            l_pts = rs.CurveCurveIntersection(long, rs.JoinCurves(bs[j]))
            lin = rs.AddLine(pt, l_pts[0][2])
            random.seed(sd*1.5)
            rand_n = random.uniform(0.35, 0.65)
            pm = rs.CurveParameter(lin, rand_n)
            
            itrs = rs.EvaluateCurve(lin, pm)
            lins = rs.SplitCurve(lin, pm)
       
            
            
            #Rotating new line
            ang = random.uniform(-15,15)
            brk = rs.RotateObject(lins[0], rs.EvaluateCurve(lins[0], pm), ang)
            
            #Create auxiliary surface to split the brep
            longbrk = rs.ExtendCurveLength(brk, 0,0, max(mx)*1000)
            longinit = rs.ExtendCurveLength(lins[1], 0,1, max(mx)*1000)
            
            sd = sd+5
            random.seed(sd)
            rnd = random.randint(0,1)
            if rnd == 0:
                new_lines = [longbrk, longinit]
                n_ln = rs.JoinCurves(new_lines)
            
            if rnd == 1:
                n_ln = long
            
            
            vc = rs.VectorCreate([0,0,0],[0,0,1])
            ln = rs.AddLine(pt, pt+vc)
            aux = rs.ExtrudeCurve(n_ln, ln)
            
            #Split the original breps
            nb = rs.SplitBrep(breps[j], aux)
            
            if nb:
                nn.append(nb)
                sd = sd/(sd+23)
            else:
                nn.append(breps[j])
                sd = sd/(sd+23)
        else:
            nn.append(breps[j])
            sd = sd/(sd+23)
            
        new_breps = []
        for ls in nn:
            if ls:
                if isinstance(ls, list):
                    for j in ls:
                        new_breps.append(j)
                else:
                    new_breps.append(ls)
            else:
                new_breps.append(breps[1])
        
        #Get the new brep edges
        nc = [rs.DuplicateEdgeCurves(bp) for bp in new_breps]
        
        #Join the sides
        ns = []
        for ls in nc:
            ns.append(rs.JoinCurves(ls))
        
    edges = [pl for ls in ns for pl in ls]
    return edges


#Function to iterate recursively
#Function to iterate recursively
def recursion(list, empty_list, iterations):
    if iterations > 0:
        list = subd(list, M, S)
        iterations -=1
        recursion(list, empty_list, iterations)
        if iterations == 0:
            empty_list.append(list)
        return empty_list


M = int(M)
I = int(I)
S= int(S)

check_curve_orientation(P)

empty = []
C = gt.list_to_tree(recursion(P, empty, I))
