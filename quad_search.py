import numpy as np
import pickle
from pdb import set_trace as keyboard
import matplotlib.pyplot as plt
from skimage.transform import ProjectiveTransform

f = open("quad_list_leaves.pickle", "rb")
quad_list_leaves, quad_list_unit, quad_list_index, quad_list_depth, uni_old, ind_file = pickle.load(f)
response = raw_input("Enter the coordinates type T-P or rho-e: ")
rootrect = [0.14, -70793, 1305.20, 600000]

b_left  = [rootrect[0],rootrect[1]]
b_right = [rootrect[2],rootrect[1]]
t_left  = [rootrect[0],rootrect[3]]
t_right = [rootrect[2],rootrect[3]]

### Transform to a square of dimensions [1024, 1024] ###
trans = ProjectiveTransform()
src = np.asarray([b_left, b_right, t_left, t_right])
dst = np.asarray([[0,0], [1024,0], [0,1024] ,[1024,1024]])
if not trans.estimate(src, dst):
    raise Exception("estimate failed")
    keyboard()
rho_list = []
Eint_list = []
points_list = []

for i,item in enumerate(quad_list_unit):
    pointp = item[0]
    point_00 = [pointp[0],pointp[1]]
    point_10 = [pointp[2],pointp[1]]
    point_01 = [pointp[0],pointp[3]]
    point_11 = [pointp[2],pointp[3]]
    points_list.append(point_00)
    points_list.append(point_10)
    points_list.append(point_01)
    points_list.append(point_11)

    if response == "T-P":
        temp_list.append(item[0][0])
        temp_list.append(item[0][2])
        pres_list.append(item[0][1])
        pres_list.append(item[0][3])
    if response == "rho-e":
        rho_list.append(item[0][0])
        rho_list.append(item[0][2])
        Eint_list.append(item[0][1])
        Eint_list.append(item[0][3])

set(rho_list)
set(Eint_list)

keyboard()
def bilinear_interpolation(x, y, points, trans):
    q00=[None]*9
    q01=[None]*9
    q10=[None]*9
    q11=[None]*9
    prop_array=[]
    ##################------here points x0, y0 are of the order BL, TL, BR, TR
    keyboard()
    #x, y = np.ndarray.tolist(trans.inverse([x, y])[0])
    print x, y
    for p in range(9):
        q00[p]  = points[0][1][p]
        x0, y0  = np.ndarray.tolist(trans.inverse(points[0][0])[0])
        q01[p]  = points[2][1][p]
        _x0, y1 = np.ndarray.tolist(trans.inverse(points[2][0])[0])
        q10[p]  = points[1][1][p]
        x1, _y0 = np.ndarray.tolist(trans.inverse(points[1][0])[0])
        q11[p]  = points[3][1][p]
        _x1, _y1= np.ndarray.tolist(trans.inverse(points[3][0])[0])

        #commented due to the roundoff error by the transforming function
        # if x0 != _x0 or x1 != _x1 or y0 != _y0 or y1 != _y1:
        #     raise ValueError('points do not form a rectangle')
        #     sys.exit()
        if not x0 <= x <= x1 or not y0 <= y <= y1:
            raise ValueError('(x, y) not within the rectangle')
            sys.exit()

        bil = (q00[p] * (x1 - x) * (y1 - y) + q10[p] * (x - x0) * (y1 - y) + q01[p] * (x1 - x) * (y - y0) + q11[p] * (x - x0) * (y - y0)) / ((x1 - x0) * (y1 - y0) + 0.0)
        prop_array.append(bil)
    return prop_array

if response == "rho-e":

	rho_x, Eint_x = [float(x) for x in raw_input("Enter the unknown temperature, Pressure [T, P] in [K, Pa] (WITHOUT BRACES): ").split(',')]

	#x, y = np.ndarray.tolist(trans([rho_x, Eint_x])[0])
	#find uniform index [idx][idy]
	x, y = rho_x, Eint_x
	print x, y
	drho = 0.99999999999
	de   = 0.99999999999
	idx  = np.floor(x/drho)
	idx  = idx.astype(int)
	idy  = np.floor((y)/de)
	idy  = idy.astype(int)

	print "indexes are, ", idx, idy
	n1 = int(uni_old[idx-1][idy-1])
	n1 = ind_file[n1+1]
	n1 = [i for i,x in enumerate(quad_list_index) if x == n1][0]

	data = quad_list_leaves[n1]
	print data
	keyboard()
	box = data[0]
	point_00 = [[box[0],box[1]], data[0]]
	point_10 = [[box[2],box[1]], data[1]]
	point_01 = [[box[0],box[3]], data[2]]
	point_11 = [[box[2],box[3]], data[3]]

	points = [point_00, point_10, point_01, point_11]

	print bilinear_interpolation(rho_x, Eint_x, points, trans)
	keyboard()

elif response == "T,P":
	f = open("quad_list_leaves.pickle", "rb")
	quad_list_leaves, quad_list_index = pickle.load(f)

	rho_x, Eint_x = [float(x) for x in raw_input("Enter the unknown Density, Internal Energy [rho, Eint] in [Kg/m3, KJ/kg] (WITHOUT BRACES): ").split(',')]

	




