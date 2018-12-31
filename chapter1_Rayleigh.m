clc; clear all

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
for i = 1:N
    for j = 1:N
        VP(i,j) = 1000;
        VS(i,j) = 570;
        Density(i,j) = 2*1000;
        lame1(i,j) = Density(i,j) * (VP(i,j)^2 - 2*VS(i,j)^2);
        lame2(i,j) = Density(i,j) * VS(i,j)^2;
    end
end

%%% U��V��P��Q��R���ٶȡ�Ӧ��ȫ�ռ��ʼ�� %%%
for i = 1:N
	for j = 1:N
		% ȫx��Ĳ���:
		U_x(i,j) = 0;
		V_x(i,j) = 0;
		P_x(i,j) = 0;
		Q_x(i,j) = 0;
		R_x(i,j) = 0;
		% ȫy��Ĳ���:
		U_y(i,j) = 0;
		V_y(i,j) = 0;
		P_y(i,j) = 0;
		Q_y(i,j) = 0;
		R_y(i,j) = 0;
		% ʱ��:
		U(i,j) = 0;
		V(i,j) = 0;
		P(i,j) = 0;
		Q(i,j) = 0;
		R(i,j) = 0;
        % 
        dx(i) = 0;
        dy(j) = 0;
	end
end

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
			if j >= N-4-layer/space_x  %  400 - 4 - 10/0.5 = 346  ÿһ��i�£�jҪ��346֮��Ż᲻����0
                %fprintf('��%d��"��ѭ��"�ѽ���j>=346�׶�,��ǰ:j=%d\n',i,j);
                %pause(1);
				dy(j) = (2*VP(i,j)/layer)*log(1/coefficientR)*((j-(N-4-layer/space_x))*space_x/layer)^4;
			end
			% �ٶȵ�������
			U_x(i,j) = ( (1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*U_x(i,j) + ...
			           (1/(1+0.5*sample_t*dx(i)))*(1/Density(i,j))*(sample_t/space_x)*...
			           ( C1*(P(i+1,j)-P((i),j))+C2*( P((i+2),j)-P((i-1),j)));  % �Ѽ��
			U_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*U_y(i,j) + ...
			           (1/(1+0.5*sample_t*dy(j)))*(1/Density(i,j))*(sample_t/space_z)* ...
			           ( C1*(R(i,j)-R(i,j-1))+C2*(R(i,(j+1))-R(i,(j-2))));
			U(i,j) = U_x(i,j) + U_y(i,j); 

			V_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*V_x(i,j) + ...
			           (1/(1+0.5*sample_t*dx(i)))*(1/Density(i,j))*(sample_t/space_x)*...
			           ( C1*(R(i,j)-R(i-1,j))+C2*(R((i+1),j)-R(i-2,j)));
			V_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*V_y(i,j) + ...
			           (1/(1+0.5*sample_t*dy(j)))*(1/Density(i,j))*(sample_t/space_z)* ...
			           ( C1*(Q(i,(j+1))-Q(i,j))+C2*( Q(i,(j+2))-Q(i,(j-1))));
			V(i,j) = V_x(i,j) + V_y(i,j);
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

    		% ������Դ��λ��:
    		V(200,4) = 100*(1-2*(pi*fm*(k-1000)*sample_t)^2)*...
    		           exp(-(pi*fm*(k-1000)*sample_t)^2);
    		% ��Դ������:
    		P_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*P_x(i,j) + ...
    		           (1/(1+0.5*sample_t*dx(i)))*(lame1(i,j)+2*lame2(i,j))*(sample_t/space_x)*...
    		           (C1*(U((i),j)-U(i-1,j))+C2*(U((i+1),j)-U((i-2),j)));
    		P_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*P_y(i,j) + ...
    		           (1/(1+0.5*sample_t*dy(j)))*lame1(i,j)*(sample_t/space_z)*...
    		           (C1*(V(i,j)-V(i,(j-1))) + C2*( V(i,(j+1))-V(i,(j-2))));
    		P(i,j) = P_x(i,j) + P_y(i,j);

    		Q_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*Q_x(i,j) + ...
    		           (1/(1+0.5*sample_t*dx(i)))*lame1(i,j)*(sample_t/space_x)*...
    		           (C1*(U((i),j)-U(i-1,j))+C2*(U((i+1),j)-U((i-2),j)));
    		Q_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*Q_y(i,j) + ...
    		           (1/(1+0.5*sample_t*dy(j)))*(lame1(i,j)+2*lame2(i,j))*(sample_t/space_z)*...
                       (C1*(V(i,j)-V(i,(j-1)))+C2*(V(i,j+1)-V(i,(j-2))));
            Q(i,j) = Q_x(i,j) + Q_y(i,j);

            R_x(i,j) = ((1-0.5*sample_t*dx(i))/(1+0.5*sample_t*dx(i)))*R_x(i,j) + ...
                       (1/(1+0.5*sample_t*dx(j)))*lame2(i,j)*(sample_t/space_x)*...
                       (C1*(V(i+1,j)-V((i),j)) + C2*( V((i+2),j)-V((i-1),j)));
            R_y(i,j) = ((1-0.5*sample_t*dy(j))/(1+0.5*sample_t*dy(j)))*R_y(i,j) + ...
                       (1/(1+0.5*sample_t*dy(j)))*lame2(i,j)*(sample_t/space_z)*...
                       (C1*(U(i,(j+1))-U(i,j)) + C2*(U(i,(j+2))-U(i,(j-1))));
            R(i,j) = R_x(i,j) + R_y(i,j);
        end
    end
    % 130msʱ�Ĳ������գ�k = t*20
    if k == 10  
        figure('name','���ٶ�V�Ĳ�������');
        colormap(gray);
        imagesc(V')
        figure('name','�ٶȷ���V_x�Ĳ�������');
        colormap(gray);
        imagesc(V_x');
        figure('name','�ٶȷ���V_y�Ĳ�������');
        colormap(gray);
        imagesc(V_y');
        break;
    end
end
toc
xlswrite('result_rayleigh.xls', V', 'sheet1')
xlswrite('result_rayleigh.xls', V_x', 'sheet2')
xlswrite('result_rayleigh.xls', V_y', 'sheet3')







       
