import pdb
import matplotlib.pyplot as plt
#from quad_plot import draw_rectangle
import pickle
import numpy as np
import NIST_reader as NIST
import quad_utilities as utility
import math
from pdb import set_trace as keyboard
from skimage.transform import ProjectiveTransform

class Node():
    ROOT = 0
    BRANCH = 1
    LEAF = 2
    minsize = 1   # Set by QuadTree
    #_______________________________________________________.
    # In the case of a root node "parent" will be None. The
    # "rect" lists the minx,minz,maxx,maxz of the rectangle
    # represented by the node.
    def __init__(self, parent, rect, rect_prop, parent_index, n, accuracy, response, trans):
        self.parent = parent
        self.children = [None,None,None,None] #initially no children
        if parent == None:
            self.depth = 0 #starting point
        else:
            self.depth = parent.depth + 1
        self.rect = rect
        self.index = (4*parent_index) + (n+1)
        self.rect_prop = rect_prop
        x0,z0,x1,z1 = rect #initial outline for the grid

        if self.parent == None:
            self.type = Node.ROOT #if there is no parent then this is the root
            self.index = 0
        elif (x1 - x0) <= Node.minsize:
            self.type = Node.LEAF # if there is a parent, see if this is the final node i.e., leaf
        else:
            self.type = Node.BRANCH #if it is not parent and not a leaf, then it is a branch
    #______________________________________________________
    # Recursively subdivides a rectangle. Division occurs 
    # ONLY if the rectangle spans a "feature of interest" or the error is large enough.
    def subdivide(self, parent_index, accuracy, response, trans):
        if self.type == Node.LEAF: #if it is a leaf you can't go any further..so return nothing
            return
        x0,z0,x1,z1 = self.rect #assign the outline coordinates to the rectangle
        dx = (x1 - x0)/2
        dz = (z1 - z0)/2

        rects = [] 
        rects_prop = []
        rects.append( (x0, z0, x0 + dx, z0 + dz) ) #just appending the list of new child rect coordinates
        rects.append( (x0, z0 + dz, x0 + dx, z1) )
        rects.append( (x0 + dx, z0 + dz, x1, z1) )
        rects.append( (x0 + dx, z0, x1, z0 + dz) )
        
        for n in range(len(rects)):
            #check if they're out of bounds
            # check = utility.check_out_of_bound_points(rects[n], response, trans)
            # if check:
            #     rects_prop.append(self.rect_prop)
            # else:
            rects_prop.append( utility.get_coolprop_TP(rects[n], response, trans))
            self.children[n] = self.getinstance(rects[n], rects_prop[n], parent_index, n, accuracy, response, trans) #assigning the class function to that child        

        # for n in range(len(rects)):
        #     #rects_prop.append( utility.get_dataNIST_I(rects[n]) )
        #     rects_prop.append( utility.get_coolprop_TP(rects[n]))
        #     span = self.spans_feature(rects[n], rects_prop[n], self.depth, accuracy) #for each child, check if it spans a feature
        #     self.children[n] = self.getinstance(rects[n], rects_prop[n], parent_index, n, accuracy) #assigning the class function to that child
        #     if span == True:
        #         #if not(self.type == Node.ROOT):
        #             keyboard()
        #             print "Subdividing further in level: ",self.depth
        #             print "Current point's index:    ", self.children[n].index
        #             self.children[n].subdivide(self.children[n].index, accuracy) # << recursion
                #else:
                    #Node.subdivide()
        for no, child in enumerate(self.children):
            span = self.spans_feature(child.rect, child.rect_prop, child.depth, accuracy, response, trans) #for each child, check if it spans a feature
            if span == True:
                    print "Subdividing further in level: ",self.depth
                    print "Current point's index:    ", self.children[n].index
                    child.subdivide(child.index, accuracy, response, trans) # << recursion



    #_______________________________________________________
    # Sub-classes must override these two methods.
    def getinstance(self,rect):
        return Node(self,rect)            

  
#===========================================================            
class QuadTree():
    maxdepth = 1 # the "depth" of the tree
    leaves = []
    allnodes = []

    #_______________________________________________________
    def __init__(self, rootnode, minrect, accuracy, response, trans):
        Node.minsize = minrect
        QuadTree.max_ref_level =  4
        rootnode.subdivide(0, accuracy, response, trans) # constructs the network of nodes
        self.prune(rootnode)
        #index = 0
        outputName = "quad_list_leaves.pickle"
        quad_list = [None]*100000
        quad_list_leaves = []
        quad_list_index = []
        quad_list_unit = []
        quad_list_depth = []
        self.traverse(rootnode, quad_list)
        # pdb.set_trace()
        #quad_list = quad_list[::-1]#reversing the original
        # for count_none in range(0, len(quad_list)):
        #     if quad_list[count_none]!= None:
        #         break
        #count_none =  #finding the first 'None' instance
        # quad_list = quad_list[::-1]
        #pdb.set_trace()
        #act = 1000000 - count_none
        #quad_list = quad_list[:act]
        #quad_list[0] = rootnode.rect

        ####----------saving data in lists---------------####
        for i,element in enumerate(QuadTree.leaves):
            prop_list = [element.rect, element.rect_prop[0], element.rect_prop[1], element.rect_prop[3], element.rect_prop[2]]
            quad_list_unit.append(prop_list)
            data = [[element.rect[0],element.rect[1]], [element.rect[2],element.rect[1]], [element.rect[0],element.rect[3]], [element.rect[2],element.rect[3]]]
            tp_data = np.ndarray.tolist(trans.inverse(data))
            tp_data = [tp_data[0][0], tp_data[0][1], tp_data[3][0], tp_data[3][1]]
            tp_prop_list = [tp_data, element.rect_prop[0], element.rect_prop[1], element.rect_prop[3], element.rect_prop[2]]
            quad_list_leaves.append(tp_prop_list)
            quad_list_index.append(element.index)
            quad_list_depth.append(element.depth)
        
        ####-----------generating uniform array index -------###
        n = 0
        uni_old = np.zeros((1,1))
        for n in range(QuadTree.maxdepth):
            uni_new = np.zeros((2**(n+1), 2**(n+1)))
            for i in range(2**n): 
                for j in range(2**n):
                    uni_new[2*i][2*j] = int(uni_old[i][j]*4 + 1)
                    uni_new[(2*i)+1][2*j] = int(uni_old[i][j]*4 + 2)
                    uni_new[(2*i)+1][(2*j)+1] = int(uni_old[i][j]*4 + 3)
                    uni_new[2*i][(2*j)+1] = int(uni_old[i][j]*4 + 4)
            uni_old = uni_new
        
        ####-----------modifying indices according to the generated quad_tree----##
        ind_file = range(0,int(np.max(uni_old))+2)

        def index_change(p_index, depth, index):
            depth += 1
            c1 = int((4*p_index) +1)
            c2 = int((4*p_index) +2)
            c3 = int((4*p_index) +3)
            c4 = int((4*p_index) +4)
            ind_file[c1] = index
            ind_file[c2] = index
            ind_file[c3] = index
            ind_file[c4] = index
            ch = [c1, c2, c3, c4]
            if depth<QuadTree.maxdepth:
                for j, child in enumerate(ch):
                    index_change(child, depth, index)

        for i, item in enumerate(quad_list_index):
            depth = quad_list_depth[i]
            global index
            index = item
            if depth<QuadTree.maxdepth:
                index_change(item, depth, index)

        print "Sucessfully created a Tree, writing it into file."
        print "The total no. of points are: ", len(quad_list_leaves)
        print "The tree maximum depth is: ", QuadTree.maxdepth
        f=open(outputName, "wb" )
        pickle.dump((quad_list_leaves, quad_list_unit, quad_list_index, quad_list_depth, uni_old, ind_file, trans, QuadTree.maxdepth), f )
        f.close()

        ans = raw_input("Would you like to test the tree ? (Y/N): ")
        while (ans == 'y' or ans == 'Y'):
            print ("#-----Testing.. Testing... !!")
            rho_x, Eint_x = [float(x) for x in raw_input("Enter the unknown Density, Internal Energy [rho, Eint] in [Kg/m3, KJ/kg] (WITHOUT BRACES): ").split(',')]
            rho_x, Eint_x = np.ndarray.tolist(trans([rho_x, Eint_x])[0])
            self.search(rootnode, rho_x, Eint_x)
            
            box = self.box
            box_prop = self.box_prop
            point_00 = [[box[0],box[1]], box_prop[0]]
            point_10 = [[box[2],box[1]], box_prop[1]]
            point_01 = [[box[0],box[3]], box_prop[2]]
            point_11 = [[box[2],box[3]], box_prop[3]]

            points = [point_00, point_10, point_01, point_11]
            keyboard()
            print self.bilinear_interpolation(rho_x, Eint_x, points)
            ans = raw_input("Would you like to test the tree ? (Y/N): ")

        # x, y = rho_x, Eint_x
        # rho_x, Eint_x = np.ndarray.tolist(trans([rho_x, Eint_x])[0])
        # drho = 0.99999999999
        # de   = 0.99999999999
        # idx  = np.floor(rho_x/drho)
        # idx  = idx.astype(int)
        # idy  = np.floor(Eint_x/de)
        # idy  = idy.astype(int)

        # print "indexes are, ", idx, idy
        # n1 = int(uni_old[idx-1][idy-1])
        # n1 = ind_file[n1]
        # #print "The box index as per overlap is, ", n1
        # keyboard()
        # n1 = [i for i,x in enumerate(quad_list_index) if x == n1][0]
        # data = quad_list_leaves[n1]
        # print data

        # keyboard()
    def search(self, node, rho_x, Eint_x):
        for num, child in enumerate(node.children):
            if self.contains(child, rho_x, Eint_x):
                if child.type == Node.LEAF:
                    print "The box index as per binary search is : ", child.index
                    self.box = child.rect
                    self.box_prop = child.rect_prop
                else:
                    self.search(child, rho_x, Eint_x)

    # A utility proc that returns True if the coordinates of
    # a point are within the bounding box of the node.
    def contains(self, node, x, z):
        x0,z0,x1,z1 = node.rect
        if x >= x0 and x <= x1 and z >= z0 and z <= z1:
            return True
        return False
    
    def bilinear_interpolation(self, x, y, points):
        q00=[None]*9
        q01=[None]*9
        q10=[None]*9
        q11=[None]*9
        prop_array=[]
        for p in range(9):
            q00[p]  = points[0][1][p]
            x0, y0  = points[0][0]
            q01[p]  = points[2][1][p]
            _x0, y1 = points[2][0]
            q10[p]  = points[1][1][p]
            x1, _y0 = points[1][0]
            q11[p]  = points[3][1][p]
            _x1, _y1= points[3][0]

            if x0 != _x0 or x1 != _x1 or y0 != _y0 or y1 != _y1:
                raise ValueError('points do not form a rectangle')
                sys.exit()
            if not x0 <= x <= x1 or not y0 <= y <= y1:
                raise ValueError('(x, y) not within the rectangle')
                sys.exit()

            bil = (q00[p] * (x1 - x) * (y1 - y) + q10[p] * (x - x0) * (y1 - y) + q01[p] * (x1 - x) * (y - y0) + q11[p] * (x - x0) * (y - y0)) / ((x1 - x0) * (y1 - y0) + 0.0)
            prop_array.append(bil)
        return prop_array


    #_______________________________________________________
    # Sets children of 'node' to None if they do not have any
    # LEAF nodes.        
    def prune(self, node):
        if node.type == Node.LEAF: #if it is of type leaf
            return 1
        leafcount = 0 #setting the leafcount = 0
        removals = []
        for child in node.children: 
            if child != None:
                leafcount += self.prune(child)
                if leafcount == 0:
                    removals.append(child)
        # for item in removals:
        #     n = node.children.index(item)
        #     node.children[n] = None        
        return leafcount
    #_______________________________________________________
    # Appends all nodes to a "generic" list, but only LEAF 
    # nodes are appended to the list of leaves.

    def traverse(self, node, quad_list):
        QuadTree.allnodes.append(node)
        if (node.type == Node.LEAF):
            QuadTree.leaves.append(node)
            if node.depth > QuadTree.maxdepth:
                QuadTree.maxdepth = node.depth
        for no, child in enumerate(node.children):
            if node.index >= 1000000:
                print "More than what we expected, please increase the list size"
                f=open(outputName, "wb" )
                pickle.dump(quad_list, f )
                f.close()
                pdb.set_trace()
            quad_list.insert(node.children[no].index, node.children[no].rect)
            #quad_list.insert(node.index, new_list) #list contains all the variables
            if node.children[no].children != [None, None, None, None]:
                self.traverse(node.children[no], quad_list) # << recursion
            else:
                node.children[no].type = Node.LEAF
                QuadTree.leaves.append(node.children[no])
                if node.children[no].depth > QuadTree.maxdepth:
                    QuadTree.maxdepth = node.children[no].depth
