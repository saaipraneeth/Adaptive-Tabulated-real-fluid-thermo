    
from Quad_tree_UI import Node, QuadTree
#from quad_plot import draw_rectangle
import sys
import random
from pdb import set_trace as keyboard
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import NIST_reader as NIST
import pickle
import quad_utilities_t as utility
import numpy as np
from CoolProp.CoolProp import PropsSI
from skimage.transform import ProjectiveTransform
  

class CNode(Node):
    #_______________________________________________________
    # Overrides the base class method.
    # Ensures Node.subdivide() uses instances of our custom 
    # class rather than instances of the base 'Node' class.
    def getinstance(self, rect, rect_prop, index, n, accuracy, response, trans, uniform_index_list):
        return CNode(self, rect, rect_prop, index, n, accuracy, response, trans, uniform_index_list)
    
    # Test if the reconstructed values are within the error limits of the acutal values, if not subdivide
    def spans_feature(self, rect, rect_prop, depth, accuracy, response, trans, uniform_index_list):
        x0,z0,x1,z1 = rect 
        print self.depth
        if self.depth>2:
            return False

        return True


class CQuadTree(QuadTree):
    #_______________________________________________________
    def __init__(self, rootnode, minrect, accuracy, response, trans, uniform_index_list):
        QuadTree.__init__(self, rootnode, minrect, accuracy, response, trans, uniform_index_list)
    

if __name__=="__main__":

    print "#######################----ADAPTIVE TABULATION PROGRAM FOR THERMODYNAMIC EQUATION OF STATE------###########################"
    print "      "
    response = "rho-e"

    #rho_min, rho_max = [float(x) for x in raw_input("Enter the range of density [rho_min, rho_max] in Kg/m3 (WITHOUT BRACES): ").split(',')]
    #e_min, e_max = [float(x) for x in raw_input("Enter the range of energies [e_min, e_max] in J/kg: (WITHOUT BRACES) ").split(',')] 
    #rootrect = [rho_min, e_min, rho_max, rho_max]
    #rootrect = [0.14,-154110, 1208, 333470]
    #rootrect = [1.4,-12, 717, 163000]
    #rootrect = [1.4,-12000, 78, 163000]
    rootrect = [0.14, -70793, 1305.20, 600000]

    b_left  = [rootrect[0],rootrect[1]]
    b_right = [rootrect[2],rootrect[1]]
    t_left  = [rootrect[0],rootrect[3]]
    t_right = [rootrect[2],rootrect[3]]

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
    print("The initial rectangular domain of rho-e is ", rootrect)

    # points = [point_00, point_10, point_01, point_11]
    
    rootrect_prop = [] #the properties at those given ranges(only on 4 boundary points)

    rootrect_prop = utility.get_coolprop_TP(rootrect, response, trans)
    uniform_index_list = []

    resolution = 1
    accuracy = 1
    rootnode = CNode(None, rootrect, rootrect_prop, 0, 0, accuracy, response, trans, uniform_index_list)
    tree = CQuadTree(rootnode, resolution, accuracy, response, trans, uniform_index_list)
    #print "Done"
    #pdb.set_trace()
    #f=open("quadtree.pickle", "wb" )
    #pickle.dump((tree), f )
    #f.close()

    print "End of program"

# initial = np.zeros((2,2))

# initial_n = 1
# initial[0][0] = 1
# initial[1][0] = 2
# initial[1][1] = 3
# initial[0][1] = 4

# n = 1
# uni_old = np.zeros((1,1))
# while n<=3:
#     uni_new = np.zeros((2**n, 2**n))
#     for i in range(2**n): 
#         for j in range(2**n):
#             uni_new[2*i][2*j] = uni_old[i][j]*4 + 1
#             uni_new[(2*i)+1][2*j] = uni_old[i][j]*4 + 2
#             uni_new[(2*i)+1][(2*j)+1] = uni_old[i][j]*4 + 3
#             uni_new[2*i][(2*j)+1] = uni_old[i][j]*4 + 4
#     uni_old = uni_new
#     n += 1
