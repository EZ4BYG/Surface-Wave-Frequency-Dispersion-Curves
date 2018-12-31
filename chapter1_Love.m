clc; clear all

tic
%%% 参数设置 %%%
N = 400;
fm = 30;
space_x = 0.5; space_z = 0.5;  % 网格步长
sample_t = 0.05/1000;  % 时间步长0.05ms
T = 250/1000;    % 总时间为250ms
K = T/sample_t;   % 外层时间循环总次数:250/0.05 = 5000次
layer = 20*space_x;  % 吸收层的厚度
coefficientR = 0.0001;  % 理想的反射系数：基本不反射

C1 = 1.125;
C2 = -0.04166667;

%%% 均匀半空间模型的初始化 %%%
for i = 1:N
    for j = 1:N
        VP(i,j) = 1000;
        VS(i,j) = 570;
        Density(i,j) = 2*1000;
        lame1(i,j) = Density(i,j) * (VP(i,j)^2 - 2*VS(i,j)^2);
        lame2(i,j) = Density(i,j) * VS(i,j)^2;
    end
end

%%% Love波:速度、应力初始化 %%%
for i = 1:N
    for j = 1:N
        % 速度参数:震源设置在vy(i,j)上即可
        vyx(i,j) = 0;
        vyz(i,j) = 0;
        vy(i,j) = 0;
        % 应力参数:
        txy(i,j) = 0;
        tzy(i,j) = 0;
        % 间距参数:
        dx(i) = 0;
        dy(j) = 0;
    end
end

%%% AEA 自由表面边界条件 %%%
% 其实就是把原始空间的前4行拿出来，用来做自由空间层！
% 所以差分计算时，i和j都是从4开始
for i = 1:4
    for j = 1:N
        Density(i,j) = 0.5*Density(i,j);
        lame1(i,j) = 0;
        lame2(i,j) = lame2(i,j);
    end
end

%%% 开始有限差分的计算 %%%
for k = 1:K  % 时间是最外层循环：相当于总的迭代次数/传播的时间
    fprintf('迭代计算开始,当前是第%d次"时间"大循环!\n',k)
    mibinbin = k;
    % AEA 自由边界设置：
    for i = 1:N
        Q(i,4) = 0;
    end

    % 差分计算：各点的速度更新
    for i = 4:(N-4)
        for j = 4:(N-4)
            % 区域边界的设置
            if i <= 4 + layer/space_x
                dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((4+layer/space_x-i)*space_x/layer)^4;
            end
            if i >= N-4-layer/space_x
                dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((i-(N-4-layer/space_x))*space_x/layer)^4;
            end
            if j >= N-4-layer/space_x 
                % 400 - 4 - 10/0.5 = 346  每一个i下，j要在346之后才会不等于0
                % fprintf('第%d次"行循环"已进入j>=346阶段,当前:j=%d\n',i,j);
                % pause(1);
                dy(j) = (2*VP(i,j)/layer)*log(1/coefficientR)*((j-(N-4-layer/space_x))*space_x/layer)^4;
            end
            % 速度迭代计算
            vyx(i,j) = (vyx(i,j)*(1-sample_t*dx(i)/2)+sample_t/Density(i,j)*((C1*(txy(i,j)-txy(i,j-1))+...
                        C2*(txy(i,j+1)-txy(i,j-2)))/space_x))/(1+sample_t*dx(i)/2);
            vyz(i,j) = (vyz(i,j)*(1-sample_t*dy(i)/2)+sample_t/Density(i,j)*((C1*(tzy(i,j)-tzy(i-1,j))+...
                        C2*(tzy(i+1,j)-tzy(i-2,j)))/space_z))/(1+sample_t*dy(i)/2);
            vy(i,j) = vyx(i,j) + vyz(i,j);
        end
    end

    % 差分计算：各点的应力更新
    for i = 4:(N-4)
        for j = 4:(N-4)
            if i <= 4 + layer/space_x
                dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((4+layer/space_x-i)*space_x/layer)^4;
            end
            if i >= N-4-layer/space_x
                dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((i-(N-4-layer/space_x))*space_x/layer)^4;
            end
            if j >= N-4-layer/space_x
                dy(j) = (2*VP(i,j)/layer)*log(1/coefficientR)*((j-(N-4-layer/space_x))*space_x/layer)^4;
            end
            % 震源位置设置:
            vy(200,4) = 100*(1-2*(pi*fm*(k-1000)*sample_t)^2)*...
                       exp(-(pi*fm*(k-1000)*sample_t)^2);
            % 应力设置:
            txy(i,j) = (lame1(i,j)*sample_t*((C1*(vy(i,j+1)-vy(i,j))+C2*(vy(i,j+2)-vy(i,j-1))/space_x))+...
                        txy(i,j)*(1-sample_t*dx(i)/2))/(1+sample_t*dx(i)/2);
            tzy(i,j) = (lame1(i,j)*sample_t*((C1*(vy(i+1,j)-vy(i,j))+C2*(vy(i+2,j)-vy(i-1,j))/space_z))+...
                        tzy(i,j)*(1-sample_t*dy(i)/2))/(1+sample_t*dy(i)/2);
        end
    end
    % 130ms时的波场快照：k = t*20
    if k == 2600 
        figure('name','Love波总速度V的波场快照');
        colormap(jet);
        imagesc(vy')
        figure('name','Love波y-x方向速度的波场快照');
        colormap(jet);
        imagesc(vyx');
        figure('name','Love波y-z方法速度的波场快照');
        colormap(jet);
        imagesc(vyz');
        break;
    end
end
toc
xlswrite('result_love.xls', vy', 'sheet1');
xlswrite('result_love.xls', vyx', 'sheet2');
xlswrite('result_love.xls', vyz', 'sheet3');