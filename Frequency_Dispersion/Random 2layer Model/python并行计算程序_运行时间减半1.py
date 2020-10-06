from multiprocessing import Process,Pipe
import numpy as np
import sympy  
import datetime

'''
说明1: 使用的是python3.5 模块最好都升级到最新版本;
说明2: 现在multiprocessing模块已经是python3.5及以上版本自带模块,无需安装!
说明3: 频率范围最好是1-240Hz;
说明4: 每个进程计算60个频率所对应的相速度/根,各进程间无重叠计算!
说明5: 进程1中有"假根"出现,即相对于已知结果来说它是错误的!但不知其错误的原由是什么,故只能通过测试剔除
'''

# 模型已知参数: 
vp1 = float(input('第一层p波速度(m/s):'))
vp2 = float(input('第二层p波速度(m/s):'))
vp = [vp1,vp2]     # 纵波速度
vs1 = float(input('第一层s波速度(m/s):'))
vs2 = float(input('第二层s波速度(m/s):'))
vs = [vs1,vs2]     # 横波速度
den1 = float(input('第一层密度(kg/m^3):'))
den2 = float(input('第二层密度(kg/m^3):'))
den = [den1,den2]
h1 = float(input('第一层厚度(m):'))
h = h1
n = 2    # 两种的介质

# 两个字符变量!!
w, Vr = sympy.symbols('w,Vr')
kk = w/Vr

rp = [sympy.sqrt(Vr*Vr/vp[0]**2-1), sympy.sqrt(Vr*Vr/vp[1]**2-1)]
rs = [sympy.sqrt(Vr*Vr/vs[0]**2-1), sympy.sqrt(Vr*Vr/vs[1]**2-1)]
r = [1-(Vr*Vr)/(2*(vs[0]**2)), 1-(Vr*Vr)/(2*(vs[1]**2))]
g = [1-r[0], 1-r[1]]
rr = [rp[0]**2, rp[1]**2]
s = [rs[0]**2, rs[1]**2]
p = [rp[0]*kk*h, rp[1]*kk*h]
q = [rs[0]*kk*h, rs[1]*kk*h]
a = [sympy.cos(p[0]), sympy.cos(p[1])]
b = [sympy.cos(q[0]), sympy.cos(q[1])]
c = [sympy.sin(p[0])/rp[0], sympy.sin(p[1])/rp[1]]
d = [sympy.sin(q[0])/rs[0], sympy.sin(q[1])/rs[1]]
l = [vs[0]**2*den[0]/(vs[0]**2*den[0]), 
     vs[1]**2*den[1]/(vs[1]**2*den[1])]

F = np.zeros([5,5])
E = np.array([[0,0,0,0,0],
              [1+rp[n-1]*rs[n-1], r[n-1]+rp[n-1]*rs[n-1], rs[n-1]*(1-r[n-1])*sympy.I,\
               rp[n-1]*(r[n-1]-1)*sympy.I, -r[n-1]**2-rp[n-1]*rs[n-1]]
             ]).T   # .T转置

for m in range(n-1,0,-1):
    M1 = np.array([[1,2,0,0,-1],
                   [r[m-1],1+r[m-1],0,0,-1],
                   [0,0,g[m-1],0,0],
                   [0,0,0,g[m-1],0],
                   [-r[m-1]**2,-2*r[m-1],0,0,1]])
    L = np.array([[a[m-1]*b[m-1],0,-a[m-1]*d[m-1],b[m-1]*c[m-1],c[m-1]*d[m-1]],
                  [0,1,0,0,0],
                  [a[m-1]*d[m-1]*s[m-1],0,a[m-1]*b[m-1],c[m-1]*d[m-1]*s[m-1],-b[m-1]*c[m-1]],
                  [-b[m-1]*c[m-1]*rr[m-1],0,c[m-1]*d[m-1]*rr[m-1],a[m-1]*b[m-1],a[m-1]*d[m-1]],
                  [c[m-1]*d[m-1]*rr[m-1]*s[m-1],0,b[m-1]*c[m-1]*rr[m-1],-a[m-1]*d[m-1]*s[m-1],a[m-1]*b[m-1]]])
    M2 = np.array([[1/l[m-1],-2,0,0,-l[m-1]],
                   [-r[m-1]/l[m-1],1+r[m-1],0,0,l[m-1]],
                   [0,0,g[m-1],0,0],
                   [0,0,0,g[m-1],0],
                   [-r[m-1]**2/l[m-1],2*r[m-1],0,0,l[m-1]]])
    F = np.dot(M1,L,M2)
    E[:,m-1] = np.dot(F,E[:,m])

fun = E[4,0]
# print(fun)
print('快速矢量传递结束!\n')

wwmin = int(input('二分求根,最小面波频率(Hz):'))
wwmax = int(input('二分求根,最大面波频率(Hz):'))

# 时间记录:
starttime = datetime.datetime.now()

# 进程1计算第一个1/4频率点
def pro1():
    # 循环体设置
    ww_min = 1
    ww_max = wwmax//4
    ww_step = 1
    Vrr_min = 0.81*min(vs)
    Vrr_max = 1.2*max(vs)
    Vrr_step = 0.55
    ww = ww_min    # 循环起始
    Vrr = Vrr_min  # 循环起始

    # 求根参数设置
    root_tmp1 = 0  # 临时1:记录根 
    root_tmp2 = 0  # 临时2:记录根
    k = 1          # 根个数计数器;若有根个数肯定从1开始!
    acc = 0.0001   # 误差/精度
    root1 = np.zeros([4,wwmax])  # 记录根的(Vrr的一系列数值)

    print('进程1:第一个1/4部分频率求根开始!')
    while ww < ww_max:
        while Vrr < Vrr_max:
            left = Vrr
            right = Vrr + Vrr_step
            gap = right - left
            if Vrr > max(vs):
                break
            # 测试: print('1:',fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right}))
            tmp = fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right})
            if ( (sympy.re(tmp) < -1*10**(-7)) and (ww <=10) ) or \
               ( (sympy.re(tmp) < -1*10**(-8)) and (ww > 10)):
                print('1:',sympy.re(tmp))
                print('进程1二分求根阶段!')
                while gap > acc:
                    center = (left+right)/2
                    if sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) < 0:
                        right = center
                    elif sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) > 0:
                        left = center
                    else:
                        left = center
                        right = center
                    gap = gap/2
                root_tmp1 = (left+right)/2
                if abs(root_tmp1-root_tmp2) > 0.55:
                    root1[k-1,ww-1] = root_tmp1
                    root_tmp2 = root1[k-1,ww-1]
                    print('进程1: 频率:%d Hz 对应的相速度/根是:%.6f m/s\n'%(ww,root1[k-1,ww-1]));
                    k = k+1
            Vrr = Vrr + Vrr_step 
        k = 1
        ww = ww + ww_step
        Vrr = Vrr_min  # 非常容易忘!
    # 求根任务结束,将根发给主进程
    print('进程1结束!')
    child1_conn.send(root1)

# 进程2计算第二个1/4频率点
def pro2():
    # 循环体设置
    ww_min = wwmax//4
    ww_max = 2*wwmax//4
    ww_step = 1
    Vrr_min = 0.81*min(vs)
    Vrr_max = 1.2*max(vs)
    Vrr_step = 0.55
    ww = ww_min    # 循环起始
    Vrr = Vrr_min  # 循环起始

    # 求根参数设置
    root_tmp1 = 0  # 临时1:记录根 
    root_tmp2 = 0  # 临时2:记录根
    k = 1          # 根个数计数器;若有根个数肯定从1开始!
    acc = 0.0001   # 误差/精度
    root2 = np.zeros([4,wwmax])  # 记录根的(Vrr的一系列数值)

    print('进程2:第二个1/4部分频率求根开始!')
    while ww < ww_max:
        while Vrr < Vrr_max:
            left = Vrr
            right = Vrr + Vrr_step
            gap = right - left
            if Vrr > max(vs):
                break
            # 测试: print('2:',fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right}))
            tmp = fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right})
            if sympy.re(tmp) < 0:
                print('进程2二分求根阶段!')
                while gap > acc:
                    center = (left+right)/2;
                    if sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) < 0:
                        right = center
                    elif sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) > 0:
                        left = center
                    else:
                        left = center
                        right = center
                    gap = gap/2
                root_tmp1 = (left+right)/2
                if abs(root_tmp1-root_tmp2) > 0.55:
                    root2[k-1,ww-1] = root_tmp1
                    root_tmp2 = root2[k-1,ww-1]
                    print('进程2: 频率:%d Hz 对应的相速度/根是:%.6f m/s\n'%(ww,root2[k-1,ww-1]));
                    k = k+1
            Vrr = Vrr + Vrr_step
        k = 1
        ww = ww + ww_step
        Vrr = Vrr_min
    # 求根任务结束,将根发给主进程
    print('进程2结束!')
    child2_conn.send(root2)

# 进程3计算第三个1/4频率点
def pro3():
    # 循环体设置
    ww_min = 2*wwmax//4
    ww_max = 3*wwmax//4
    ww_step = 1
    Vrr_min = 0.81*min(vs)
    Vrr_max = 1.2*max(vs)
    Vrr_step = 0.55
    ww = ww_min    # 循环起始
    Vrr = Vrr_min  # 循环起始

    # 求根参数设置:
    root_tmp1 = 0  # 临时1:记录根 
    root_tmp2 = 0  # 临时2:记录根
    k = 1          # 根个数计数器;若有根个数肯定从1开始!
    acc = 0.0001   # 误差/精度
    root3 = np.zeros([4,wwmax])  # 记录根的(Vrr的一系列数值)

    print('进程3:第三个1/4部分频率求根开始!')
    while ww < ww_max:
        while Vrr < Vrr_max:
            left = Vrr
            right = Vrr + Vrr_step
            gap = right - left
            if Vrr > max(vs):
                break
            # print('3:',fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right}))
            tmp = fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right})
            if sympy.re(tmp) < 0:
                print('进程3二分求根阶段!')
                while gap > acc:
                    center = (left+right)/2;
                    if sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) < 0:
                        right = center
                    elif sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) > 0:
                        left = center
                    else:
                        left = center
                        right = center
                    gap = gap/2
                root_tmp1 = (left+right)/2
                if abs(root_tmp1-root_tmp2) > 0.55:
                    root3[k-1,ww-1] = root_tmp1
                    root_tmp2 = root3[k-1,ww-1]
                    print('进程3: 频率:%d Hz 对应的相速度/根是:%.6f m/s\n'%(ww,root3[k-1,ww-1]));
                    k = k+1
            Vrr = Vrr + Vrr_step
        k = 1
        ww = ww + ww_step
        Vrr = Vrr_min
    # 求根任务结束,将根发给主进程
    print('进程3结束!')
    child3_conn.send(root3)

# 进程4计算第四个1/4频率点
def pro4():
    # 循环体设置
    ww_min = 3*wwmax//4
    ww_max = 4*wwmax//4
    ww_step = 1
    Vrr_min = 0.81*min(vs)
    Vrr_max = 1.2*max(vs)
    Vrr_step = 0.55
    ww = ww_min    # 循环起始
    Vrr = Vrr_min  # 循环起始

    # 求根参数设置
    root_tmp1 = 0  # 临时1:记录根 
    root_tmp2 = 0  # 临时2:记录根
    k = 1          # 根个数计数器;若有根个数肯定从1开始!
    acc = 0.0001   # 误差/精度
    root4 = np.zeros([4,wwmax])  # 记录根的(Vrr的一系列数值)

    print('进程4:第四个1/4部分频率求根开始!')
    while ww < ww_max: 
        while Vrr < Vrr_max:
            left = Vrr
            right = Vrr + Vrr_step
            gap = right - left
            if Vrr > max(vs):
                break
            # 测试: print('4:',fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right}))
            tmp = fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:right})
            if sympy.re(tmp) < 0:
                print('进程4二分求根阶段!')
                while gap > acc:
                    center = (left+right)/2;
                    if sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) < 0:
                        right = center
                    elif sympy.re(fun.evalf(subs={w:ww,Vr:left})*fun.evalf(subs={w:ww,Vr:center})) > 0:
                        left = center
                    else:
                        left = center
                        right = center
                    gap = gap/2
                root_tmp1 = (left+right)/2
                if abs(root_tmp1-root_tmp2) > 0.55:
                    root4[k-1,ww-1] = root_tmp1
                    root_tmp2 = root4[k-1,ww-1]
                    print('进程4: 频率:%d Hz 对应的相速度/根是:%.6f m/s\n'%(ww,root4[k-1,ww-1]));
                    k = k+1
            Vrr = Vrr + Vrr_step
        k = 1
        ww = ww + ww_step
        Vrr = Vrr_min
    # 求根任务结束,将根发给主进程
    print('进程4结束!')
    child4_conn.send(root4)


# 回到主进程:
# 管道设置: 在进程设置之前
child1_conn, parent1_conn = Pipe()
child2_conn, parent2_conn = Pipe()
child3_conn, parent3_conn = Pipe()
child4_conn, parent4_conn = Pipe()

# 进程设置:
p1 = Process(target = pro1)
p2 = Process(target = pro2)
p3 = Process(target = pro3)
p4 = Process(target = pro4)

p1.start()
p2.start()
p3.start()
p4.start()

p1.join()
p2.join()
p3.join()
p4.join()

# 进程间通信:管道
root1 = parent1_conn.recv()
root2 = parent2_conn.recv()
root3 = parent3_conn.recv()
root4 = parent4_conn.recv()
print('所有进程二分求根结束!')

# 时间记录:
endtime = datetime.datetime.now()
print((starttime-endtime).seconds)

# 把各进程计算的部分进行合并(提取有效段+矩阵水平拼接)
# 有个完整的root就可以在matlab里画图
root_end1 = np.hstack(root1[:,0:wwmax//4],root2[:,wwmax//4:2*wwmax//4])
root_end2 = np.hstack(root3[:,2*wwmax//4:3*wwmax//4],root4[:,3*wwmax//4:wwmax])
root = hstack(root_end1,root_end2)
