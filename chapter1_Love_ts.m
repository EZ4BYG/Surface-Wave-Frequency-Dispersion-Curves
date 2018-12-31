function []= chapter1_Love_ts()

clc; clear all

snapshot_time = input('������һ����������ʱ��(��λ:ms):')

tic
%%% �������� %%%
N = 400;
fm = 30;
space_x = 0.5; space_z = 0.5;  % ���񲽳�
sample_t = 0.05/1000;  % ʱ�䲽��0.05ms
T = 250/1000;    % ��ʱ��Ϊ250ms
K = T/sample_t;   % ���ʱ��ѭ���ܴ���:250/0.05 = 5000��
layer = 20*space_x;  % ���ղ�ĺ��
coefficientR = 0.0001;  % ����ķ���ϵ��������������

C1 = 1.125;
C2 = -0.04166667;

%%% ���Ȱ�ռ�ģ�͵ĳ�ʼ�� %%%
% Ԥ��������:
VP = zeros(400,400)+1000;
VS = zeros(400,400)+570;
Density = zeros(400,400)+2000;
lame1 = zeros(400,400);
lame2 = zeros(400,400);
for i = 1:N
    j = 1:N;  % ʸ��������
    lame1(i,j) = Density(i,j) .* (VP(i,j).^2 - 2*VS(i,j).^2);
    lame2(i,j) = Density(i,j) .* VS(i,j).^2;
end

%%% Love��:�ٶȡ�Ӧ����ʼ�� %%%
% �ٶȷ���:
vyx = zeros(400,400);
vyz = zeros(400,400);
vy = zeros(400,400);
% Ӧ������:
txy = zeros(400,400);
tzy = zeros(400,400);
% ���:
dx = zeros(1,400);
dy = zeros(1,400);

%%% AEA ���ɱ���߽����� %%%
% ��ʵ���ǰ�ԭʼ�ռ��ǰ4���ó��������������ɿռ�㣡
% ���Բ�ּ���ʱ��i��j���Ǵ�4��ʼ
for i = 1:4
    for j = 1:N
        Density(i,j) = 0.5*Density(i,j);
        lame1(i,j) = 0;
        lame2(i,j) = lame2(i,j);
    end
end

%%% ��ʼ���޲�ֵļ��� %%%
for k = 1:K  % ʱ���������ѭ�����൱���ܵĵ�������/������ʱ��
    fprintf('�������㿪ʼ,��ǰ�ǵ�%d��"ʱ��"��ѭ��!\n',k)
    mibinbin = k;
    % AEA ���ɱ߽����ã�
    for i = 1:N
        Q(i,4) = 0;
    end

    % ��ּ��㣺������ٶȸ���
    for i = 4:(N-4)
        for j = 4:(N-4)
            % ����߽������
            if i <= 4 + layer/space_x
                dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((4+layer/space_x-i)*space_x/layer)^4;
            end
            if i >= N-4-layer/space_x
                dx(i) = (2*VP(i,j)/layer)*log(1/coefficientR)*((i-(N-4-layer/space_x))*space_x/layer)^4;
            end
            if j >= N-4-layer/space_x 
                % 400 - 4 - 10/0.5 = 346  ÿһ��i�£�jҪ��346֮��Ż᲻����0
                % fprintf('��%d��"��ѭ��"�ѽ���j>=346�׶�,��ǰ:j=%d\n',i,j);
                % pause(1);
                dy(j) = (2*VP(i,j)/layer)*log(1/coefficientR)*((j-(N-4-layer/space_x))*space_x/layer)^4;
            end
            % �ٶȵ�������
            vyx(i,j) = (vyx(i,j)*(1-sample_t*dx(i)/2)+sample_t/Density(i,j)*((C1*(txy(i,j)-txy(i,j-1))+...
                        C2*(txy(i,j+1)-txy(i,j-2)))/space_x))/(1+sample_t*dx(i)/2);
            vyz(i,j) = (vyz(i,j)*(1-sample_t*dy(i)/2)+sample_t/Density(i,j)*((C1*(tzy(i,j)-tzy(i-1,j))+...
                        C2*(tzy(i+1,j)-tzy(i-2,j)))/space_z))/(1+sample_t*dy(i)/2);
            vy(i,j) = vyx(i,j) + vyz(i,j);
        end
    end

    % ��ּ��㣺�����Ӧ������
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
            % ��Դλ������:
            vy(200,4) = 100*(1-2*(pi*fm*(k-1000)*sample_t)^2)*...
                       exp(-(pi*fm*(k-1000)*sample_t)^2);
            % Ӧ������:
            txy(i,j) = (lame1(i,j)*sample_t*((C1*(vy(i,j+1)-vy(i,j))+C2*(vy(i,j+2)-vy(i,j-1))/space_x))+...
                        txy(i,j)*(1-sample_t*dx(i)/2))/(1+sample_t*dx(i)/2);
            tzy(i,j) = (lame1(i,j)*sample_t*((C1*(vy(i+1,j)-vy(i,j))+C2*(vy(i+2,j)-vy(i-1,j))/space_z))+...
                        tzy(i,j)*(1-sample_t*dy(i)/2))/(1+sample_t*dy(i)/2);
        end
    end
    % snapshotʱ�Ĳ������գ�k = snapshot_time*20
    if k == snapshot_time*20 
        figure('name','Love�����ٶ�V�Ĳ�������');
        colormap(jet);
        imagesc(vy')
        figure('name','Love��y-x�����ٶȵĲ�������');
        colormap(jet);
        imagesc(vyx');
        figure('name','Love��y-z�����ٶȵĲ�������');
        colormap(jet);
        imagesc(vyz');
        break;
    end
end
toc
xlswrite('result_love.xls', vy', 'sheet1');
xlswrite('result_love.xls', vyx', 'sheet2');
xlswrite('result_love.xls', vyz', 'sheet3');


