	
from quadtree import Node, QuadTree
#from quad_plot import draw_rectangle
import sys
import random
from pdb import set_trace as keyboard
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import NIST_reader as NIST
import pickle
import quad_utilities as utility
import numpy as np
from CoolProp.CoolProp import PropsSI
from skimage.transform import ProjectiveTransform

class CNode(Node):
    #_______________________________________________________
    # Overrides the base class method.
    # Ensures Node.subdivide() uses instances of our custom 
    # class rather than instances of the base 'Node' class.
    def getinstance(self, rect, rect_prop, index, n, accuracy, response, trans):
        return CNode(self, rect, rect_prop, index, n, accuracy, response, trans)
    
    # Test if the reconstructed values are within the error limits of the acutal values, if not subdivide
    def spans_feature(self, rect, rect_prop, depth, accuracy, response, trans):
        x0,z0,x1,z1 = rect 

        # dataNIST=NIST.readNIST(isoType = "isotherm", fluid = 'O2', T=x_mid, P=z_mid/1.0E6, tmin=x_mid, tmax=x_mid, pmin = z_mid/1.0E6, pmax = z_mid/1.0E6, N=1)
        
        # dataNIST = np.ndarray.tolist(dataNIST)
        # dataNIST = [item for sublist in dataNIST for item in sublist]
        # del dataNIST[:2]
        #dataNIST = utility.get_coolprop_TPS(x_mid, z_mid, response)

        point_00 = [[rect[0],rect[1]], rect_prop[0]]
        point_10 = [[rect[2],rect[1]], rect_prop[1]]
        point_01 = [[rect[0],rect[3]], rect_prop[2]]
        point_11 = [[rect[2],rect[3]], rect_prop[3]]

        points = [point_00, point_10, point_01, point_11]

        #reconstructing property values at num x num points inside each box
        num = 7.0
        del_x = (x1-x0)/(num-1.0)
        del_z = (z1-z0)/(num-1.0)
        x_l = np.linspace(x0+del_x,x1-del_x,num-2.0)
        z_l = np.linspace(z0+del_z,z1-del_z,num-2.0)
        xv = np.ndarray.flatten(np.meshgrid(x_l,z_l)[0])
        zv = np.ndarray.flatten(np.meshgrid(x_l,z_l)[1])

        int_points = zip(xv, zv)

        mid_prop = [None]*25
        dataNIST = [None]*25
        for i,int_point in enumerate(int_points):
            x_c = int_point[0]
            z_c = int_point[1]
            mid_prop[i] = self.bilinear_interpolation(x_c, z_c, points, trans)
            dataNIST[i] = utility.get_coolprop_TPS(x_c, z_c, response, trans)

        glob_error = [None]*25
        for n in range(25):
            lc_error = [None]*7
            for e in range(7):
                lc_error[e] = abs(mid_prop[n][e] - dataNIST[n][e])/dataNIST[n][e]
            glob_error[n] = max(lc_error)

        #if all(item<(accuracy/100.0) for item in error) or depth >= 13:
        if (max(glob_error)<(accuracy/100.0)) or depth >= 13:
            return False
        # if depth >= 50:
        #     print error_rho, "This is the error in density"
        #     return False
        # if response=='t':
        #     return True
        # else:
        #     return False
        # response = raw_input("Response")
        # if response == 'Y':
        return True


    def bilinear_interpolation(self, x, y, points, trans):
        q00=[None]*9
        q01=[None]*9
        q10=[None]*9
        q11=[None]*9
        prop_array=[]
        ##################------here points x0, y0 are of the order BL, TL, BR, TR
        for p in range(9):
            q00[p]  = points[0][1][p]
            x0, y0  = points[0][0]
            q01[p]  = points[2][1][p]
            _x0, y1 = points[2][0]
            q10[p]  = points[1][1][p]
            x1, _y0 = points[1][0]
            q11[p]  = points[3][1][p]
            _x1, _y1= points[3][0]
        # x, y = np.ndarray.tolist(trans.inverse([x, y])[0])
        # for p in range(9):
        #     q00[p]  = points[0][1][p]
        #     x0, y0  = np.ndarray.tolist(trans.inverse(points[0][0])[0])
        #     q01[p]  = points[2][1][p]
        #     _x0, y1 = np.ndarray.tolist(trans.inverse(points[2][0])[0])
        #     q10[p]  = points[1][1][p]
        #     x1, _y0 = np.ndarray.tolist(trans.inverse(points[1][0])[0])
        #     q11[p]  = points[3][1][p]
        #     _x1, _y1= np.ndarray.tolist(trans.inverse(points[3][0])[0])
            if x0 != _x0 or x1 != _x1 or y0 != _y0 or y1 != _y1:
                raise ValueError('points do not form a rectangle')
                sys.exit()
            if not x0 <= x <= x1 or not y0 <= y <= y1:
                raise ValueError('(x, y) not within the rectangle')
                sys.exit()

            bil = (q00[p] * (x1 - x) * (y1 - y) + q10[p] * (x - x0) * (y1 - y) + q01[p] * (x1 - x) * (y - y0) + q11[p] * (x - x0) * (y - y0)) / ((x1 - x0) * (y1 - y0) + 0.0)
            prop_array.append(bil)
        return prop_array


class CQuadTree(QuadTree):
    #_______________________________________________________
    def __init__(self, rootnode, minrect, accuracy, response, trans):
        QuadTree.__init__(self, rootnode, minrect, accuracy, response, trans)
    

if __name__=="__main__":

    print "#######################----ADAPTIVE TABULATION PROGRAM FOR THERMODYNAMIC EQUATION OF STATE------###########################"
    print "      "
    response = raw_input("Enter the table index variables, T-P or rho-e:  ")

    if response == "T-P":
        #T_min, T_max = [float(x) for x in raw_input("Enter the range of temperatures [T_min, T_max] in K (WITHOUT BRACES): ").split(',')]
        #P_min, P_max = [float(x) for x in raw_input("Enter the range of temperatures [P_min, P_max] in Pa: (WITHOUT BRACES) ").split(',')]
        #rootrect = [T_min, P_min, T_max, P_max]
        rootrect = [200, 0.02E5, 600, 10.0E6]

        b_left  = [rootrect[0],rootrect[1]]
        b_right = [rootrect[2],rootrect[1]]
        t_left  = [rootrect[0],rootrect[3]]
        t_right = [rootrect[2],rootrect[3]]
        raw_rootrect = rootrect
        ### Transform to a square of dimensions [1024, 1024] ###
        trans = ProjectiveTransform()
        src = np.asarray([b_left, b_right, t_left, t_right])
        dst = np.asarray([[0,0], [1023,0], [0,1023] ,[1024,1024]])
        rootrect = [0, 0, 1024, 1024]
        if not trans.estimate(src, dst):
            raise Exception("estimate failed")
            keyboard()

        point_00 = [rootrect[0],rootrect[1]]
        point_10 = [rootrect[2],rootrect[1]]
        point_01 = [rootrect[0],rootrect[3]]
        point_11 = [rootrect[2],rootrect[3]]
        print("The initial rectangular domain of PT is", raw_rootrect)
        print("The transformed domain of PT is", rootrect)
    
    # points = [point_00, point_10, point_01, point_11]
    elif response == "rho-e":
        #rho_min, rho_max = [float(x) for x in raw_input("Enter the range of density [rho_min, rho_max] in Kg/m3 (WITHOUT BRACES): ").split(',')]
        #e_min, e_max = [float(x) for x in raw_input("Enter the range of energies [e_min, e_max] in J/kg: (WITHOUT BRACES) ").split(',')] 
        #rootrect = [rho_min, e_min, rho_max, rho_max]

        #rootrect = [0.14, -70793, 1305.20, 300000]
        rootrect = [0.14, 50000, 130.20, 300000]
        b_left  = [rootrect[0],rootrect[1]]
        b_right = [rootrect[2],rootrect[1]]
        t_left  = [rootrect[0],rootrect[3]]
        t_right = [rootrect[2],rootrect[3]]
        raw_rootrect = rootrect
        ### Transform to a square of dimensions [1024, 1024] ###
        trans = ProjectiveTransform()
        src = np.asarray([b_left, b_right, t_left, t_right])
        dst = np.asarray([[0,0], [1024,0], [0,1024] ,[1024,1024]])
        rootrect = [0, 0, 1024, 1024]
        if not trans.estimate(src, dst):
            raise Exception("estimate failed")
            keyboard()

        point_00 = [rootrect[0],rootrect[1]]
        point_10 = [rootrect[2],rootrect[1]]
        point_01 = [rootrect[0],rootrect[3]]
        point_11 = [rootrect[2],rootrect[3]]
        print("The initial rectangular domain of PT is", raw_rootrect)
        print("The transformed domain of PT is", rootrect)

    else:
        print "The variables not recongnized !!"
        sys.exit()
    # points = [point_00, point_10, point_01, point_11]
    
    rootrect_prop = [] #the properties at those given ranges(only on 4 boundary points)

    rootrect_prop = utility.get_coolprop_TP(rootrect, response, trans)


    resolution = 1
    accuracy = float(raw_input("Enter the required accuracy in (%) "))
    rootnode = CNode(None, rootrect, rootrect_prop, 0, 0, accuracy, response, trans)
    tree = CQuadTree(rootnode, resolution, accuracy, response, trans)
    #print "Done"
    #pdb.set_trace()
    #f=open("quadtree.pickle", "wb" )
    #pickle.dump((tree), f )
    #f.close()


    print "End of program"
        

