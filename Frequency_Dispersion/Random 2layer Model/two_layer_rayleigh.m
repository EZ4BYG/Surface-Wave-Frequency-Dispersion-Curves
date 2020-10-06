%%% 两层模型瑞雷波频散曲线 %%%

clc; clear;

% 模型已知参数: 自定义
% 代码总体差距不大，只需在下面数值判断的地方稍微调整即可满足自定义两层模型参数
% 测试还不完全：希望输入的数值是符合客观规律的！即P波速度大的层S波速度也是大的！
%              否则阶数/方程根变化我还不能确定！
fprintf('任意两层模型参数自定义设置:\n');
vp1 = input('第一层p波速度(m/s):');
vp2 = input('第二层p波速度(m/s):');
vp = [vp1,vp2]; 
vs1 = input('第一层s波速度(m/s):');
vs2 = input('第二层s波速度(m/s):');
vs = [vs1,vs2];  % 横波速度
den1 = input('第一层密度(km/m^3):');
den2 = input('第二层密度(km/m^3):');
den = [den1,den2];
h1 = input('第一层厚度(m):');
h = [h1];
% Vr = 200;  % 初始相速度
n = 2;  % 两种的介质

% 相速度和波数
syms w Vr;
kk = w/Vr;

% 快速矢量传递算法中间参数:计算后全是矩阵！
rp = sqrt(Vr^2./vp.^2-1);  
rs = sqrt(Vr^2./vs.^2-1);  
r = 1-Vr^2./(2*(vs.^2));   
g = 1-r;      
rr = rp.^2;   
s = rs.^2;
p = rp*kk.*h;
q = rs*kk.*h;
a = cos(p);
b = cos(q);
c = sin(p)./rp;
d = sin(q)./rs;
l = vs.^2.*den./(vs.^2.*den);

F = zeros(5,5);
%E = zeros(5,n);
% 以第二个元素作为起始
E = [0,0,0,0,0;
     1+rp(n)*rs(n), r(n)+rp(n)*rs(n), rs(n)*(1-r(n))*i, rp(n)*(r(n)-1)*i, -r(n)^2-rp(n)*rs(n)]';  
for m = n-1:-1:1
    M1 = [1,2,0,0,-1;
          r(m),1+r(m),0,0,-1;
          0,0,g(m),0,0;
          0,0,0,g(m),0;
          -r(m)^2,-2*r(m),0,0,1];
    L = [a(m)*b(m),0,-a(m)*d(m),b(m)*c(m),c(m)*d(m);
         0,1,0,0,0;
         a(m)*d(m)*s(m),0,a(m)*b(m),c(m)*d(m)*s(m),-b(m)*c(m);
         -b(m)*c(m)*rr(m),0,c(m)*d(m)*rr(m),a(m)*b(m),a(m)*d(m);
         c(m)*d(m)*rr(m)*s(m),0,b(m)*c(m)*rr(m),-a(m)*d(m)*s(m),a(m)*b(m)];
    M2 = [1/l(m),-2,0,0,-l(m);
          -r(m)/l(m),1+r(m),0,0,l(m);
          0,0,g(m),0,0;
          0,0,0,g(m),0;
          -r(m)^2/l(m),2*r(m),0,0,l(m)];
    F = M1*L*M2;
    E(:,m) = F*E(:,m+1);  % 第m列列向量：对"所有行"的第m列进行赋值:右边应该是5个数值的列向量
end
fun = E(5,1);
fprintf('快速矢量传递结束!\n');

% 这里是可以通用！
% fun = E(5,1)
% ww = 10; Vrr = 189.1;
% fun = eval(subs(fun,[w,Vr],[ww,Vrr]))

% 二分法求fun的根
wwmin = input('二分求根,最小面波频率(Hz):');
wwmax = input('二分求根,最大面波频率(Hz):');
root_tmp1 = 0;  % 临时1:记录根 
root_tmp2 = 0;  % 临时2:记录根
k = 1;          % 根个数计数器;若有根个数肯定从1开始!
acc = 0.0001;   % 误差/精度
root = zeros(4,wwmax);  % 记录根的(Vrr的一系列数值)

fprintf('二分求根开始! (频率步长默认1 Hz;速度步长默认0.44 m/s)\n');
tic % 后面我想做并行计算,这里记录单进程/线程用时
for ww = wwmin:wwmax
    % 测试发现规律1: 29Hz之后有双根(频率精确)!
    % 测试发现规律2: 150Hz之后有三根(准确三根Hz暂不定)！
    % 测试发现规律3: 200Hz之后有四根(准确四根Hz暂不定)！
    % 测试发现规律4: 250Hz之后有五根(准确五根Hz暂不定)!
    % 测试发现规律5: 240Hz四根！
    % 测试后建议: 1 - 240 Hz (两层模型预计运算时长都在6000秒左右) 
    fprintf('当前求根频率为:%d Hz\n',ww)
    % 注意: 速度步长最影响运行时间!
    for Vrr = 0.81*min(vs): 0.44 : 1.2*max(vs)
        left = Vrr;
        right = Vrr+0.44;
        gap = right - left;  % 在外面；步长横为0.44
        % 事实上不到1.2倍就可以提前结束了(跳出当前Vrr循环)，后面不会有根了
        if Vrr > max(vs)
		    break;
        end
        
        % 测试语句:
        % if eval(subs(fun,[w,Vr],[ww,left]))*eval(subs(fun,[w,Vr],[ww,right])) < 0
        %     fprintf('通过')
        %     pause(3)
        % end
        
        % 若进入到这里说明有根！理论上每个w都有根！
        if eval(subs(fun,[w,Vr],[ww,left]))*eval(subs(fun,[w,Vr],[ww,right])) < 0  % 有根
            % fprintf('已经进入当前频率的求根阶段! 进入速度为:%.5f\n',Vrr);
        	while gap>acc
                center = (left+right)/2;
		        if eval(subs(fun,[w,Vr],[ww,left]))*eval(subs(fun,[w,Vr],[ww,center]))<0
		            right = center;
		        elseif eval(subs(fun,[w,Vr],[ww,left]))*eval(subs(fun,[w,Vr],[ww,center]))>0
		            left = center;
		        else
		            left = center;
		            right = center;
		        end
		        gap = gap/2;
            end
            % 循环已达要求，根的数值可以输出了:
            root_tmp1 = (left+right)/2;
            % 如果确实是多根！列元素的行值加1！
            % 根之间大于步长，才是新根！否则不记录/不打印结果(小步长产生的假根)！
            if abs(root_tmp1-root_tmp2) > 0.44      	
                root(k,ww) = root_tmp1;
                root_tmp2 = root(k,ww);
		        fprintf('频率:%d Hz  对应的相速度/根是:%.6f m/s\n', ww, root(k,ww));
                k = k+1;
            end
		    % pause(1);
		end
		% 根判别段落结束
    end
    k = 1;
end
fprintf('二分求根完毕!\n');
toc

% 散点绘图
figure()
axis([0 100 150 450])
% 频率w做横坐标  速度Vr做纵坐标
for ww = wwmin:wwmax
    if root(1,ww) ~= 0
    	basic = plot(ww, root(1,ww), '.b');
    	hold on;
    end
end
hold on
for ww = wwmin:wwmax
    if root(2,ww) ~= 0
    	first = plot(ww, root(2,ww), '.r');
    	hold on;
    end
end
hold on
for ww = wwmin:wwmax
    if root(3,ww) ~= 0
    	second = plot(ww, root(3,ww),'.g');
    	hold on;
    end
end
hold on
for ww = wwmin:wwmax
    if root(4,ww) ~= 0
    	third = plot(ww, root(4,ww), '.m');
    	hold on;
    end
end
legend([basic,first,second,third],'基阶','一阶高阶','二阶高阶','三阶高阶')

% 将参数写入文件保存:
xlswrite('2layer-rayleigh.xlsx', root, 'sheet1');




