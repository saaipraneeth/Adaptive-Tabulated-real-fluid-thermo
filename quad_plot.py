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
import quad_utilities as utility
import itertools
from CoolProp.CoolProp import PropsSI, PhaseSI
import numpy as np
import matplotlib.cm as cm
from skimage.transform import ProjectiveTransform
import unicodedata

fluid = 'Oxygen'
response = raw_input("Enter the coordinates type T-P or rho-e: ")
f = open("quad_list_leaves.pickle", "rb")
quad_list_leaves, quad_list_unit, quad_list_index, quad_list_depth, uni_old, ind_file = pickle.load(f)
f = open("quad_list_rect.pickle", "rb")
quad_list_rect, quad_list_uindex = pickle.load(f)
#quad_list = [x for x in quad_list if x is not None]
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
        rho_list_u.append(item[0][0])
        rho_list_u.append(item[0][2])
        Eint_list_u.append(item[0][1])
        Eint_list_u.append(item[0][3])
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
        rho_list.append(item[0][0])
        rho_list.append(item[0][2])
        Eint_list.append(item[0][1])
        Eint_list.append(item[0][3])


trans = ProjectiveTransform()

if response == "T-P":
    b_left = [200.0,0.02E6]
    b_right = [600,0.02E6]
    t_left = [200,10.0E6]
    t_right = [600,10.0E6]
    src = np.asarray([b_left, t_left, t_right, b_right])
    dst = np.asarray([[0,0], [0, 1024], [1024, 1024], [1024,0]])
    if not trans.estimate(src, dst): raise Exception("estimate failed")
elif response == "rho-e":
    b_left = [0.14, -70793]
    b_right = [1305.20, -70793]
    t_left = [0.14, 600000]
    t_right = [1305.20, 600000]    
    src = np.asarray([b_left, t_left, t_right, b_right])
    dst = np.asarray([[0,0], [0, 1024], [1024, 1024], [1024,0]])
    if not trans.estimate(src, dst): raise Exception("estimate failed")


points_list.sort()
points = list(points_list for points_list,_ in itertools.groupby(points_list))

for i, item in enumerate(points):

    #rhoe_list.append([PropsSI('DMASS', 'P', item[1],'T', item[0], fluid), PropsSI('UMASS', 'P', item[1],'T', item[0], fluid)/1000.0])
    #rho_list.append(PropsSI('DMASS', 'P', item[1],'T', item[0], fluid))
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

if response == "T-P":
    currentAxis.set_xlim([min(temp_list),max(temp_list)])
    currentAxis.set_ylim([min(pres_list),max(pres_list)])
    min_color = min(x for x in color_list if x is not None)
    max_color = np.max(color_list)
    my_cmap = cm.get_cmap('jet')
    norm = matplotlib.colors.Normalize(min_color, max_color)
elif response == "rho-e":
    currentAxis.set_xlim([min(rho_list),max(rho_list)])
    currentAxis.set_ylim([min(Eint_list),max(Eint_list)])
    min_color = min(x for x in color_list if x is not None)
    max_color = np.max(color_list)
    my_cmap = cm.get_cmap('jet')
    norm = matplotlib.colors.Normalize(min_color, max_color)


# def draw_rectangle(currentAxis, node, depth):
#     """Recursively plot a visualization of the quad tree region"""
#     if depth is None or depth == 0:
#         #rect = plt.Rectangle(self.mins, *self.sizes, zorder=2, ec='#000000', fc='none')
#         #fig.add_patch(rect)
#         x1 = node.rect[0]
#         x2 = node.rect[2]
#         y1 = node.rect[1]
#         y2 = node.rect[3]
#         size_x = (x2-x1)
#         size_y = (y2-y1)
#         #currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1))
#         #figureName = figsuffix + ".png"
#         #plt.savefig(figureName,dpi=75)
#         currentAxis.plt.scatter([x1,x2,x1,x2],[y1,y2,y2,y1])
#     if depth is None or depth > 0:
#         for child in node.children:
#             draw_rectangle(currentAxis, node, depth- 1)

for i, item in enumerate(quad_list_leaves):
    if item==None:
        continue
    x1 = item[0][0]
    x2 = item[0][2]
    y1 = item[0][1]-6.0
    y2 = item[0][3]-6.0
    size_x = (x2 - x1)
    size_y = (y2 - y1)
    mid_point = [(x1+x2)/2, (y1+y2)/2]
    if response == "rho-e":
        mid_point_uc = trans(mid_point)[0]
        check = utility.check_out_of_bound(mid_point_uc[0], mid_point_uc[1], response, trans)
        if not(check):
            color = PropsSI('T', 'UMASS', mid_point[1],'DMASS', mid_point[0], fluid)
            color_i = my_cmap(norm(color))
            currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))
    elif response == "T-P":
        mid_point_uc = trans(mid_point)[0]
        # check = utility.check_out_of_bound(mid_point_uc[0], mid_point_uc[1], response, trans)
        # if not(check):
        color = PropsSI('DMASS', 'P', mid_point[1],'T', mid_point[0], fluid)
        color_i = my_cmap(norm(color))
        currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))
        plt.text(mid_point[0], mid_point[1], quad_list_index[i], color=color_i )
cmmapable = cm.ScalarMappable(norm, my_cmap)
cmmapable.set_array(range(int(min_color), int(max_color)))

if response == "T-P":
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    cbar = plt.colorbar(cmmapable)
    #cbar.set_label(r'\textbf{Density} (\rho)',fontsize=10)
    cbar.set_label(r'$\textbf{Density} [Kg/m^3]$', fontsize = 13)
    figureName = figsuffix + ".png"
    print "Saving figure: " + figureName

    plt.xlabel(r'\textbf{Temperature} [K]',fontsize=12)
    plt.ylabel(r'\textbf{Pressure} [Pa]',fontsize=12)
    plt.title(r"\textit{Adaptive look up table for Oxygen with 0.1\% accuracy}",
              fontsize=12, color='black')
elif response == "rho-e":
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    cbar = plt.colorbar(cmmapable)
    #cbar.set_label(r'\textbf{Density} (\rho)',fontsize=10)
    cbar.set_label(r'$\textbf{Pressure} [Pa]$', fontsize = 13)
    figureName = figsuffix + ".png"
    print "Saving figure: " + figureName

    plt.xlabel(r'\textbf{Density} [Kg/m^3]',fontsize=12)
    plt.ylabel(r'\textbf{Internal Energy} [J/Kg]',fontsize=12)
    plt.title(r"\textit{Adaptive look up table for Oxygen with 1\% accuracy}",
              fontsize=12, color='black')

# plt.savefig(figureName, dpi=1000)
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
            color_i = my_cmap(norm(color))
            currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))
    elif response == "T-P":
        # check = utility.check_out_of_bound(mid_point[0], mid_point[1], response, trans)
        # if not(check):
        color = PropsSI('DMASS', 'P', mid_point_regc[1],'T', mid_point_regc[0], fluid)
        color_i = my_cmap(norm(color))
        currentAxis.add_patch(Rectangle((x1, y1), size_x, size_y, fill=None, alpha=1, color=color_i))
        plt.text(mid_point[0], mid_point[1], quad_list_index[i], color=color_i )
cmmapable = cm.ScalarMappable(norm, my_cmap)
cmmapable.set_array(range(int(min_color), int(max_color)))
cbar = plt.colorbar(cmmapable)
plt.savefig("new.png", dpi = 600)
plt.show()