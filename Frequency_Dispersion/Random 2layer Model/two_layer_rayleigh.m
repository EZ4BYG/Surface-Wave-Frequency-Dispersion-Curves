%%% ����ģ�����ײ�Ƶɢ���� %%%

clc; clear;

% ģ����֪����: �Զ���
% ���������಻��ֻ����������ֵ�жϵĵط���΢�������������Զ�������ģ�Ͳ���
% ���Ի�����ȫ��ϣ���������ֵ�Ƿ��Ͽ͹۹��ɵģ���P���ٶȴ�Ĳ�S���ٶ�Ҳ�Ǵ�ģ�
%              �������/���̸��仯�һ�����ȷ����
fprintf('��������ģ�Ͳ����Զ�������:\n');
vp1 = input('��һ��p���ٶ�(m/s):');
vp2 = input('�ڶ���p���ٶ�(m/s):');
vp = [vp1,vp2]; 
vs1 = input('��һ��s���ٶ�(m/s):');
vs2 = input('�ڶ���s���ٶ�(m/s):');
vs = [vs1,vs2];  % �Შ�ٶ�
den1 = input('��һ���ܶ�(km/m^3):');
den2 = input('�ڶ����ܶ�(km/m^3):');
den = [den1,den2];
h1 = input('��һ����(m):');
h = [h1];
% Vr = 200;  % ��ʼ���ٶ�
n = 2;  % ���ֵĽ���

% ���ٶȺͲ���
syms w Vr;
kk = w/Vr;

% ����ʸ�������㷨�м����:�����ȫ�Ǿ���
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
% �Եڶ���Ԫ����Ϊ��ʼ
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
    E(:,m) = F*E(:,m+1);  % ��m������������"������"�ĵ�m�н��и�ֵ:�ұ�Ӧ����5����ֵ��������
end
fun = E(5,1);
fprintf('����ʸ�����ݽ���!\n');

% �����ǿ���ͨ�ã�
% fun = E(5,1)
% ww = 10; Vrr = 189.1;
% fun = eval(subs(fun,[w,Vr],[ww,Vrr]))

% ���ַ���fun�ĸ�
wwmin = input('�������,��С�沨Ƶ��(Hz):');
wwmax = input('�������,����沨Ƶ��(Hz):');
root_tmp1 = 0;  % ��ʱ1:��¼�� 
root_tmp2 = 0;  % ��ʱ2:��¼��
k = 1;          % ������������;���и������϶���1��ʼ!
acc = 0.0001;   % ���/����
root = zeros(4,wwmax);  % ��¼����(Vrr��һϵ����ֵ)

fprintf('���������ʼ! (Ƶ�ʲ���Ĭ��1 Hz;�ٶȲ���Ĭ��0.44 m/s)\n');
tic % �������������м���,�����¼������/�߳���ʱ
for ww = wwmin:wwmax
    % ���Է��ֹ���1: 29Hz֮����˫��(Ƶ�ʾ�ȷ)!
    % ���Է��ֹ���2: 150Hz֮��������(׼ȷ����Hz�ݲ���)��
    % ���Է��ֹ���3: 200Hz֮�����ĸ�(׼ȷ�ĸ�Hz�ݲ���)��
    % ���Է��ֹ���4: 250Hz֮�������(׼ȷ���Hz�ݲ���)!
    % ���Է��ֹ���5: 240Hz�ĸ���
    % ���Ժ���: 1 - 240 Hz (����ģ��Ԥ������ʱ������6000������) 
    fprintf('��ǰ���Ƶ��Ϊ:%d Hz\n',ww)
    % ע��: �ٶȲ�����Ӱ������ʱ��!
    for Vrr = 0.81*min(vs): 0.44 : 1.2*max(vs)
        left = Vrr;
        right = Vrr+0.44;
        gap = right - left;  % �����棻������Ϊ0.44
        % ��ʵ�ϲ���1.2���Ϳ�����ǰ������(������ǰVrrѭ��)�����治���и���
        if Vrr > max(vs)
		    break;
        end
        
        % �������:
        % if eval(subs(fun,[w,Vr],[ww,left]))*eval(subs(fun,[w,Vr],[ww,right])) < 0
        %     fprintf('ͨ��')
        %     pause(3)
        % end
        
        % �����뵽����˵���и���������ÿ��w���и���
        if eval(subs(fun,[w,Vr],[ww,left]))*eval(subs(fun,[w,Vr],[ww,right])) < 0  % �и�
            % fprintf('�Ѿ����뵱ǰƵ�ʵ�����׶�! �����ٶ�Ϊ:%.5f\n',Vrr);
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
            % ѭ���Ѵ�Ҫ�󣬸�����ֵ���������:
            root_tmp1 = (left+right)/2;
            % ���ȷʵ�Ƕ������Ԫ�ص���ֵ��1��
            % ��֮����ڲ����������¸������򲻼�¼/����ӡ���(С���������ļٸ�)��
            if abs(root_tmp1-root_tmp2) > 0.44      	
                root(k,ww) = root_tmp1;
                root_tmp2 = root(k,ww);
		        fprintf('Ƶ��:%d Hz  ��Ӧ�����ٶ�/����:%.6f m/s\n', ww, root(k,ww));
                k = k+1;
            end
		    % pause(1);
		end
		% ���б�������
    end
    k = 1;
end
fprintf('����������!\n');
toc

% ɢ���ͼ
figure()
axis([0 100 150 450])
% Ƶ��w��������  �ٶ�Vr��������
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
legend([basic,first,second,third],'����','һ�׸߽�','���׸߽�','���׸߽�')

% ������д���ļ�����:
xlswrite('2layer-rayleigh.xlsx', root, 'sheet1');




