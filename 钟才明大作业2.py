# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 01:33:32 2021

@author: 86183
"""
'''
1) 和以前一样，根据 X 和 Y 样本，确定 T"。将此值称为 T'。
2) 从 N(Vx, Cx) 生成样本 X*。 然后从 N(Vx*, Cx*) 生成 Y*sample。 从 X* 计算检验统计量和 Y* 采样并称其值为 T'。
3) 再重复Step2 W一2次，给出T, 。 .●，Tw_ ;
4) 令 R 为 T' 在 W T' 值中的排列顺序。 如果 R≤αW，则拒绝 X 在 ax 水平上的正态性的原假设。
'''
import numpy as np
import matplotlib.pyplot as plt


def mst_prim(sides,points,maps):
    mst = []
    if points <=0 or sides<points-1:
        return mst
    seleted_node = [0]						#初始化入围节点0，为起始节点，可任选
    candidate_node = [i for i in range(1, points)]	#生成候选节点(去除已入围节点)
    while len(candidate_node) > 0:			#当候选节点不为空时循环
        begin, end, minweight = 0, 0, 9999 	#初始化边的头节点、尾节点、最小权重
        for i in seleted_node:				#遍历头节点在入围节点、尾节点在候选节点的边
            for j in candidate_node:		
                if maps[i][j] < minweight:	#若当前边的权重<最小权重，则
                    minweight = maps[i][j]	#更新最小权重值
                    begin = i 					#更新边的头节点
                    end = j 					#更新边的尾节点
        mst.append([begin, end, minweight])	#将头节点、尾节点、最小权重添加到最小生成树中
        seleted_node.append(end)			#将当前尾节点添加到入围节点中
        candidate_node.remove(end)			#从候选节点中移除当前尾节点，然后继续循环
    return mst
                                            #本prim算法参照了CSDN博主「夜空下的凝视」的文章.
def my_maps(list):
    #每个点到其他点的距离作为权值
    maps=[]
    for i in Y_new:
        t=[]
        for j in Y_new:
            #将与自身的权值记为极大值，我设置为9999
            if i==j:
                t.append(9999)
            else:
                #计算两点之间的距离，即为相应的权值
                t.append(((i[0]-j[0])**2+(i[1]-j[1])**2)**0.5)
        maps.append(t)
    return maps
def my_points(maps):
    #结点数=maps数组的长度
    points=len(maps)
    return points

def my_sides(points):
    #边数=结点数-1
    sides=points-1
    return sides
def mst_T(mst):
    #求不同数据集之间的边
    #在maps列表中前一百个值为新生成的Y数据集的点的距离，后一百个则是X数据集的点的距离
    #由maps生成的mst中的首尾节点则符合这个道理
    #在MST列表中，如果头结点<=100并且尾节点>100说明是跨越两个数据集的边，或者头节点>=100,尾节点<100同理
    t=[]
    for i in range(len(mst)):
        if (mst[i][0]<=100 and mst[i][1]>100) or (mst[i][0]>100 and mst[i][1]<=100):
            x=mst[i]
            t.append(x)
            T=len(t)
    return T
def my_C(mst,points):
    #遍历结点,并存储
    points_sum=[0 for i in range(points)]
   # print(res,len(res))
   #寻找公共边,点每找到一个边值+1
    for i in mst:
        points_sum[i[0]]+=1
        points_sum[i[1]]+=1
    #print(res)
    Sum=0
    #找到每个共享公共节点的边对的数量，执行 d_i*(d_i-1)
    for i in range(len(points_sum)):
        Sum+=points_sum[i]*(points_sum[i]-1)
    #得到C值，C=Sum
    Sum=Sum/2      
    return Sum 
    

def my_T_new(T):
    #执行论文中的计算公式
    M=100
    N=100
    L=200
    T=T
    #print(T)
    C=my_C(mst,points)
    #print(C)
    E_T=2*M*N/L
    #print(E_T)  
    a=2*M*N/(L*(L-1))
    #print(a)
    b=(2*M*N-L)/L
    #print(b) 
    c=(C-L+2)/((L-2)*(L-3))
    #print(c)
    d=L*(L-1)-4*M*N+2
    #print(d) 
    Var_T_C=a*(b+c*d)
    #print(Var_T_C)
    T_new=(T-E_T)/(Var_T_C**0.5) 
    #print(T_new)
    return T_new

data_1 = 'd:\python代码\data1.txt'
data_2 = 'd:\python代码\data2.txt'
#从文件读取数据集的数据
for file in ['d:\python代码\数据集1.txt','d:\python代码\数据集2.txt']:
    with open(file=file,mode='r',encoding='utf-8') as f:
        data = f.readlines()  #txt中所有字符串读入data
        #样本1化为列表取出
        data_1_x = []
        data_1_y = []
        for line in data:
            line = line.split()        #将单个数据分隔开存好
            data_1_x.append(eval(line[0]))
            data_1_y.append(eval(line[1]))
    #两组数据集均值
    data_mean_1_x=(sum(data_1_x)/len(data_1_x))
    data_mean_1_y=(sum(data_1_y)/len(data_1_y))
    
    
    #生成新的数据集 Y,即前面的X*
    x,y = np.array(data_1_x),np.array(data_1_y)
    mean_x = [data_mean_1_x,data_mean_1_y]
    cov_x = np.cov(x,y)
    Y = np.random.multivariate_normal(mean_x, cov_x, 100)
    #生成新的数据集Y_new，即前面步骤的Y*
    Y =Y.tolist()
    Y_X=[Y[i][0] for i in range(len(Y))] 
    Y_Y=[Y[i][1] for i in range(len(Y))] 
    mean_Y=[np.mean(Y_X),np.mean(Y_Y)]
    cov_Y=np.cov((Y_X,Y_Y))#具体计算
    Y_new = np.random.multivariate_normal(mean_Y, cov_Y, 100)
    Y_new=Y_new.tolist()
        #合并X*和Y*
    Y_new.extend(Y)
    maps = my_maps(Y_new)     #权值
    points = my_points(maps)    #节点
    sides = my_sides(points)    #边
    mst=mst_prim(points,sides,maps)     #最小生成数
    T=mst_T(mst)        #X_TO_Y的连接计数
    #T'
    T_0=my_T_new(T) 
    
    #进入w循环
    w=60
    T1=[]
    for i in range(w-1):
        #生成新的数据集 Y,即前面的X*
        Y = np.random.multivariate_normal(mean_x, cov_x, 100)
        #生成新的数据集Y_new，即前面步骤的Y*
        Y =Y.tolist()
        Y_X=[Y[i][0] for i in range(len(Y))] 
        Y_Y=[Y[i][1] for i in range(len(Y))] 
        mean_Y=[np.mean(Y_X),np.mean(Y_Y)]
        cov_Y=np.cov((Y_X,Y_Y))#具体计算
        Y_new = np.random.multivariate_normal(mean_Y, cov_Y, 100)
        Y_new=Y_new.tolist()
        #合并X*和Y*
        Y_new.extend(Y)
        maps = my_maps(Y_new)     #权值
        points = my_points(maps)    #节点
        sides = my_sides(points)    #边
        mst=mst_prim(points,sides,maps)     #最小生成数
        T=mst_T(mst)        #X_TO_Y的连接计数
        #T_i即步骤中的T''
        T_i=my_T_new(T) 
        T_w=T_i
        T1.append(T_w)
    T1.append(T_0)
    #将T'放入T''中进行排列
    T1.sort()
    #print(T1)
    index=T1.index(T_0)
    if index <= 0.05*w:
        print('T_0={},index={},假设检验不通过'.format(T_0,index))
    else:
        print('T_0={},index={},假设检验通过'.format(T_0,index))
    
                
        
   