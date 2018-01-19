import numpy as np
from pdb import set_trace as keyboard
import pickle

# f = open("quad_list_leaves.pickle", "rb")
# quad_list_leaves, quad_list_unit, quad_list_index, quad_list_depth, uni_old = pickle.load(f)

# n = 0
# uni_old = np.zeros((1,1))
# for n in range(10):
# 	uni_new = np.zeros((2**(n+1), 2**(n+1)))
# 	for i in range(2**n): 
# 		for j in range(2**n):
# 			uni_new[2*i][2*j] = int(uni_old[i][j]*4 + 1)
# 			uni_new[(2*i)+1][2*j] = int(uni_old[i][j]*4 + 2)
# 			uni_new[(2*i)+1][(2*j)+1] = int(uni_old[i][j]*4 + 3)
# 			uni_new[2*i][(2*j)+1] = int(uni_old[i][j]*4 + 4)
# 	uni_old = uni_new
# keyboard()
# ind_file = range(0,int(np.max(uni_old))+2)
# def index_change(p_index, depth, index):
# 	depth += 1
# 	c1 = int((4*p_index) +1)
# 	c2 = int((4*p_index) +2)
# 	c3 = int((4*p_index) +3)
# 	c4 = int((4*p_index) +4)
# 	ind_file[c1] = index
# 	ind_file[c2] = index
# 	ind_file[c3] = index
# 	ind_file[c4] = index
# 	temp_ind_file.append(c1)
# 	temp_ind_file.append(c2)
# 	temp_ind_file.append(c3)
# 	temp_ind_file.append(c4)
# 	ch = [c1, c2, c3, c4]
# 	if depth<10:
# 		for j, child in enumerate(ch):
# 			index_change(child, depth, index)

# for i, item in enumerate(quad_list_index):
# 	depth = quad_list_depth[i]
# 	global index
# 	index = item
# 	global temp_ind_file
# 	temp_ind_file = []
# 	if depth<10:
# 		index_change(item, depth, index)
	#keyboard()


	# c2 = int((4*p_index) +2)
	# c3 = int((4*p_index) +3)
	# c4 = int((4*p_index) +4)
	# ind_file[c2] = c2
	# ind_file[c3] = c3
	# ind_file[c4] = c4
