from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.pylab import *
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import pickle
from pdb import set_trace as keyboard
import quad_utilities_t as utility
import itertools
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.cm as cm
from skimage.transform import ProjectiveTransform

fluid = 'Oxygen'
response = raw_input("Enter the coordinates type T-P or rho-e: ")
f = open("quad_list_leaves.pickle", "rb")
quad_list_leaves, quad_list_unit, quad_list_index, quad_list_depth, uni_old, ind_file, trans, max_depth= pickle.load(f)
f = open("quad_list_rect.pickle", "rb")
quad_list_rect, quad_list_uindex = pickle.load(f)
figsuffix = "/home/administrator/Desktop/Tabulated_EoS/plot"
temp_list = []
pres_list = []
points_list = []
temp_ulist = []
pres_ulist = []
points_u_list = []
rhoe_list = []
rho_list = []
Eint_list = []
tp_list = []
indexes = []
color_list = []
rho_list_u = []
Eint_list_u = []


for i,item in enumerate(quad_list_leaves):
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
for i,item in enumerate(quad_list_unit):
    pointp = item[0]
    point_00 = [pointp[0],pointp[1]]
    point_10 = [pointp[2],pointp[1]]
    point_01 = [pointp[0],pointp[3]]
    point_11 = [pointp[2],pointp[3]]
    points_u_list.append(point_00)
    points_u_list.append(point_10)
    points_u_list.append(point_01)
    points_u_list.append(point_11)

    if response == "T-P":
        temp_ulist.append(item[0][0])
        temp_ulist.append(item[0][2])
        pres_ulist.append(item[0][1])
        pres_ulist.append(item[0][3])
    if response == "rho-e":
        rho_list_u.append(item[0][0])
        rho_list_u.append(item[0][2])
        Eint_list_u.append(item[0][1])
        Eint_list_u.append(item[0][3])


#t, p = [float(x) for x in raw_input("Enter the unknown Temperature, Pressure [T, P] in [K, Pa] (WITHOUT BRACES): ").split(',')]


def bilinear_interpolation(x, y, points):
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

def get_tableprops(rho, Eint_mass, key):
    drho = 1024/(2**max_depth)
    de   = 1024/(2**max_depth)
    rho, Eint_mass = np.ndarray.tolist(trans([rho, Eint_mass])[0])
    idx  = np.floor((rho-0)/drho)
    idx  = idx.astype(int)
    idy  = np.floor((Eint_mass-0)/de)
    idy  = idy.astype(int)

    n1 = int(uni_old[idy][idx])
    n1 = int(ind_file[n1])
    n1 = [i for i,x in enumerate(quad_list_index) if x == n1][0]
    data = quad_list_unit[n1]

    rect = data[0] # rho-e rectangle / box in which the point lies
    rect_prop = data[1:] # data stored in those four corner points

    point_00 = [[rect[0],rect[1]], rect_prop[0]]
    point_10 = [[rect[2],rect[1]], rect_prop[1]]
    point_01 = [[rect[0],rect[3]], rect_prop[3]]
    point_11 = [[rect[2],rect[3]], rect_prop[2]]

    points = [point_00, point_10, point_01, point_11]
    proparray = bilinear_interpolation(rho, Eint_mass, points)
    if key == 'all':
        return proparray
    elif key == 'T':
        return proparray[0]
    elif key == 'P':
        return proparray[1]
    elif key == 'a':
        return proparray[6]
    else:
        return None

def getT_from_rhop(rho_target, P_target):
    '''
    Univariate Newton Solver, given any one of the rhoe, outputs the other variables
    '''
    import scipy
    def func(eint):
      return get_tableprops(rho_target, eint, 'P') - P_target
    #scipy.optimize.newton
    E_target = scipy.optimize.newton(func, 194000, fprime=None, args=(), tol=1.48e-08, maxiter=50, fprime2=None)
    T_target = get_tableprops(rho_target, E_target, 'T')
    return T_target,E_target


ans = raw_input("Would you like to test the tree ? (Y/N): ")
while (ans == 'y' or ans == 'Y'):
    print ("#-----Testing.. Testing... !!")
    rho_x, Eint_x = [float(x) for x in raw_input("Enter the unknown Density, Internal Energy [rho, Eint] in [Kg/m3, KJ/kg] (WITHOUT BRACES): ").split(',')]

    proparray = get_tableprops(rho_x, Eint_x, 'all')
    print "Here's the properties list, ", proparray
    ans = raw_input("Would you like to test the tree ? (Y/N): ")

keyboard()
points_list.sort()
points = list(points_list for points_list,_ in itertools.groupby(points_list))
for i, item in enumerate(points):

    if response == "T-P":
        #item = np.ndarray.tolist(trans([item[0], item[1]])[0])
        color_list.append(PropsSI('DMASS', 'P', item[1],'T', item[0], fluid))
    elif response == "rho-e":
        #check if they're out of bounds
        #item = np.ndarray.tolist(trans([item[0], item[1]])[0])
        check = utility.check_out_of_bound(item[0], item[1], response, trans)
        if check:
            color_list.append(None)
        else:
            color_list.append(PropsSI('T', 'UMASS', item[1],'DMASS', item[0], fluid))
    #tp_list.append([PropsSI('P', 'UMASS', item[1],'DMASS', item[0], fluid), PropsSI('T', 'UMASS', item[1],'DMASS', item[0], fluid)])
    #currentAxis.scatter(PropsSI('T', 'UMASS', item[1],'DMASS', item[0], fluid), PropsSI('P', 'UMASS', item[1],'DMASS', item[0], fluid))
    #currentAxis.scatter(item[0], item[1])


fig = plt.figure()
currentAxis = plt.gca()


# if response == "T-P":
#     currentAxis.set_xlim([min(temp_list),max(temp_list)])
#     currentAxis.set_ylim([min(pres_list),max(pres_list)])
#     min_color = min(x for x in color_list if x is not None)
#     max_color = np.max(color_list)
#     my_cmap = cm.get_cmap('jet')
#     norm = matplotlib.colors.Normalize(min_color, max_color)
# elif response == "rho-e":
#     currentAxis.set_xlim([min(rho_list),max(rho_list)])
#     currentAxis.set_ylim([min(Eint_list),max(Eint_list)])
#     min_color = min(x for x in color_list if x is not None)
#     max_color = np.max(color_list)
#     my_cmap = cm.get_cmap('jet')
#     norm = matplotlib.colors.Normalize(min_color, max_color)



# for i, item in enumerate(quad_list_leaves):
#     if item==None:
#         continue
#     x1 = item[0][0]
#     x2 = item[0][2]
#     y1 = item[0][1]
#     y2 = item[0][3]
#     size_x = (x2 - x1)
#     size_y = (y2 - y1)
#     mid_point = [(x1+x2)/2, (y1+y2)/2]
#     if response == "rho-e":
#         mid_point_uc = trans(mid_point)[0]
#         check = utility.check_out_of_bound(mid_point_uc[0], mid_point_uc[1], response, trans)
#         if not(check):
#             color = PropsSI('T', 'UMASS', mid_point[1],'DMASS', mid_point[0], fluid)
#             color_i = my_cmap(norm(color))
#             currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))
#     elif response == "T-P":
#         mid_point_uc = trans(mid_point)[0]
#         # check = utility.check_out_of_bound(mid_point_uc[0], mid_point_uc[1], response, trans)
#         # if not(check):
#         color = PropsSI('DMASS', 'P', mid_point[1],'T', mid_point[0], fluid)
#         color_i = my_cmap(norm(color))
#         currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))

# cmmapable = cm.ScalarMappable(norm, my_cmap)
# cmmapable.set_array(range(int(min_color), int(max_color)))

# if response == "T-P":
#     plt.rc('text', usetex=True)
#     plt.rc('font', family='serif')
#     cbar = plt.colorbar(cmmapable)
#     #cbar.set_label(r'\textbf{Density} (\rho)',fontsize=10)
#     cbar.set_label(r'$\textbf{Density} [Kg/m^3]$', fontsize = 13)
#     figureName = figsuffix + ".png"
#     print "Saving figure: " + figureName

#     plt.xlabel(r'\textbf{Temperature} [K]',fontsize=12)
#     plt.ylabel(r'\textbf{Pressure} [Pa]',fontsize=12)
#     plt.title(r"\textit{Adaptive look up table for Oxygen with 0.1\% accuracy}",
#               fontsize=12, color='black')
# elif response == "rho-e":
#     plt.rc('text', usetex=True)
#     plt.rc('font', family='serif')
#     cbar = plt.colorbar(cmmapable)
#     #cbar.set_label(r'\textbf{Density} (\rho)',fontsize=10)
#     cbar.set_label(r'$\textbf{Pressure} [Pa]$', fontsize = 13)
#     figureName = figsuffix + ".png"
#     print "Saving figure: " + figureName

#     plt.xlabel(r'\textbf{Density} [Kg/m^3]',fontsize=12)
#     plt.ylabel(r'\textbf{Internal Energy} [J/Kg]',fontsize=12)
#     plt.title(r"\textit{Adaptive look up table for Oxygen with 1\% accuracy}",
#               fontsize=12, color='black')

#plt.savefig(figureName, dpi=1000)
# plt.show()

plt.clf()
fig = plt.figure()
currentAxis = plt.gca()
currentAxis.set_xlim([0, 1024])
currentAxis.set_ylim([0, 1024])

for i, item in enumerate(quad_list_unit):
    if item==None:
        continue

    x1 = item[0][0]
    x2 = item[0][2]
    y1 = item[0][1]
    y2 = item[0][3]
    size_x = (x2 - x1)
    size_y = (y2 - y1)

    mid_point = [(x1+x2)/2, (y1+y2)/2]
    mid_point_regc = np.ndarray.tolist(trans.inverse(mid_point)[0])
    if response == "rho-e":
        check = utility.check_out_of_bound(mid_point[0], mid_point[1], response, trans)
        if not(check):
            color = PropsSI('T', 'UMASS', mid_point_regc[1],'DMASS', mid_point_regc[0], fluid)
            #color_i = my_cmap(norm(color))
            currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1))
            plt.text(mid_point[0], mid_point[1], quad_list_index[i])
    elif response == "T-P":
        # check = utility.check_out_of_bound(mid_point[0], mid_point[1], response, trans)
        # if not(check):
        color = PropsSI('DMASS', 'P', mid_point_regc[1],'T', mid_point_regc[0], fluid)
        color_i = my_cmap(norm(color))
        currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))
#cmmapable = cm.ScalarMappable(norm, my_cmap)
#cmmapable.set_array(range(int(min_color), int(max_color)))
#cbar = plt.colorbar(cmmapable)
plt.savefig("new.png", dpi = 600)