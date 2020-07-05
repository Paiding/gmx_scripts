# coding: utf-8

import os
#import numpy as np
import shutil

path = "box.gro"
save_path = "cluster.gro"
repr_atom = 21 #用第21个原子代表整个分子的位置
mole_num = 1000 #可以读完行数算出来，但没必要
atom_num = 42 #每个分子含多少原子
MinPts = 4 
eps = 2.0

def readFile(file):
	text_ = []
	structrue = {}
	with open(file,'r') as f:
		text_ = f.readlines()
	#line_num = text_[1]
	for i in xrange(0,mole_num):#覆盖原子坐标部分，倒数第二行是line_num+1
		n = 42 * i + 2 + repr_atom
		if n < 10000:
			list_ = text_[n].strip().split()[3:6]
		else:
			list_ = text_[n].strip().split()[2:5]
		structrue[i] = list_
	return structrue

def neibourCheck(k, structrue, num, eps):
	#返回某个点的邻域内所有点的编号
	x0 = float(structrue[k][0]) - eps
	x1 = float(structrue[k][0]) + eps
	y0 = float(structrue[k][1]) - eps
	y1 = float(structrue[k][1]) + eps
	z0 = float(structrue[k][2]) - eps
	z1 = float(structrue[k][2]) + eps
	neib_list = []
	for i in range(0, num - 1):
		x = float(structrue[i][0])
		y = float(structrue[i][1])
		z = float(structrue[i][2])
		if ((x0 < x) & (x < x1)):
			if ((y0 < y) & (y < y1)):
				if ((z0 < z) & (z < z1)):
					neib_list.append(i)

	return neib_list

def unProcessedNum(state, nlist):
	u_num = 0
	for i in nlist:
		if (state[i] == 0):
			u_num = u_num + 1
	return u_num

def DBScan(structrue, num, MinPts, eps):
	cluster = {}
	c_num = 0 #给聚类计数 
	clist = [] #一个列表，表示各个聚类中点数的多少
	state = [0] * num#0未处理，1边界点，2噪声点，3核心点
	
	for i in range(1,num):
		if (state[i] > 0):
			continue
		else:
			nlist = neibourCheck(i, structrue, num, eps)
			n_num = len(nlist)
			if (n_num > MinPts):
				state[i] = 3#标记为核心点
				print("New cluster found")
				r_num = unProcessedNum(state, nlist) #该中心对应的边界点还有多少没Check过？初始化
				while r_num > 1:#这里按说最好是0，不知道为啥总是死循环，待修复
					for j in range(0, n_num - 1):
						k = nlist[j] #被检查的边界点对应的编号
						if (state[k] == 0):
							k_nlist = neibourCheck(k, structrue, num, eps)#点k
							nlist = nlist + k_nlist#聚类吸收
							tmp = set(nlist)#去重
							nlist = list(tmp)
							state[k] = 1#标记为边界点
							n_num = len(nlist)
					r_num = unProcessedNum(state, nlist) #更新剩余未处理点的数目，是否跳出循环
				clist.append(len(nlist))
				header = [i]
				cluster[c_num] = (header + nlist)
				c_num += 1
			elif (len(nlist) == 0):
				state[i] = 2#标记为噪声点
			else:
				state[i] = 1#标记为边界点
	
	return cluster, clist, c_num

st = readFile(path)

cl , clist, c_num = DBScan(st, mole_num, MinPts, eps)

with open(save_path, mode = 'w') as file:
	for i in range(0, c_num - 1):
		file.write(str(cl) + '\n')

print("complete")

