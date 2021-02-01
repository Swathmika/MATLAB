clear all;
close all;
tic;
filename='newdata.xlsx';
h = 8.77; %heat constant in W/m2K
%Te= 294:1:300; %Environment/Surrounding Temperature 
Te = 298;
R=0.001:0.001:0.005; % radius of tumour 
thetha=163.70; %angle at the center of the breast
r=0.03:0.01:0.07;  % radius of Breast
arclen=r*thetha*(pi/180);  %formula for arclength calcluation
a = -arclen/2:0.001:arclen/2 ; %defining the points and range on the breast surface
% a=0;
Q = 0.02: 0.02:0.1; %variable heat intensity 
d = 0.01:0.002:0.018 ; %fixed depth 
n=1;
Datax={'Te','Q','d','R','r','T','a'}
xlswrite(filename,Datax,'Range','A1')
%for i=1:1:length(Te)
    for j=1:1:length(r)
        for k=1:1:length(Q)
            for l=1:1:length(R)
                for m=1:1:length(d)
                 for i= 1:1:length(a)  
                    T = Te + (Q(k) /(4*h*(((d(m)+R(l))^2)+(a(i)^2))));%tempertaure at the surface calculation
                    T=T- 273.15; 
   % Data(n,1)=Te-273;
                    Data(n,2)=Q(k);
                    Data(n,3)=d(m);
                    Data(n,4)=R(l);
                    Data(n,5)=r(j);
                    Data(n,6)=T;
                    temp(n,1)=T;
                    Data(n,7)=a(i);
                    n=n+1;
                 end
                end 
            end
        end
    end
 

figure
plot(a,temp(1:86)); % 2D plot
xlabel('Location(m)'); %it's the 'a' parameter that the effected breast surface distance 
ylabel('Temperature(°C)'); % temperature at the surface
