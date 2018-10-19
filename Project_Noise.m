%PROJECT%
clear variables
close all
clc

%VARIABLES%
C = 0.9;
tau1 = 0.4;
tau2 = 0.15;
snr = 15;


%%
%PLOTTING THE DISPLACEMENT EQUATION%

n = 300;
i = 1;
for t = linspace(-4,4,n);
    if t<0;
        x(i) = abs(cos(pi/2*(t-tau1))).*exp(-(C*sin(2*pi*(t-tau1)))).*tan(pi./(1+exp(-(t-tau1))));      %Equation when t is negative
    elseif t > 0;
        x(i) = abs(cos(pi/2*(t+tau2))).*exp(-(C*sin(2*pi*(t+tau2)))).*tan(pi./(1+exp(-(t+tau2))));      %Equation when t is positive
    else
        x(i) = 0;       %Equation when t is zero
    end
    i = i+1;
end
t = linspace(-4,4,n);
x_noise = add_awgn_noise(x,snr);
plot(t,x,'Color',[.1,.1,.34]);
hold on;
plot(t,x_noise,'Color',[.2,.6,.15]);
hold on;

%%
%GLOBAL FIT%

PolDegree = 36;         %Degree of polynomial
n = PolDegree + 1;      %Number of points in fit

t = linspace(-4,4,n);
t = t';

%For-loop discretizes equation and stores values in vector x
for i = 1:length(t);
    if t(i) < 0;
        x_2 = abs(cos(pi/2*(t(i)-tau1))).*exp(-(C*sin(2*pi*(t(i)-tau1)))).*tan(pi./(1+exp(-(t(i)-tau1))));     %Equation when t is negative
        x_noise2(i) = add_awgn_noise(x_2,snr);      %Adding noise
    elseif t(i) > 0;
        x_2 = abs(cos(pi/2*(t(i)+tau2))).*exp(-(C*sin(2*pi*(t(i)+tau2)))).*tan(pi./(1+exp(-(t(i)+tau2))));     %Equation when t is positive
        x_noise2(i) = add_awgn_noise(x_2,snr);      %Adding noise
    else
        x_2 = 0;       %Equation when t is zero
        x_noise2(i) = add_awgn_noise(x_2,snr);      %Adding noise
    end
end

%Approach to solve polynomial fit: a = (T^t*T^-1)^(-1)*T^t*x
T_matrix = [];
    for k = 1:length(t)
        T_matrix = [T_matrix t.^(k-1)];
    end
% T_matrix = fliplr(vander(t));       %Creting T matrix
% a = inv(transpose(T_matrix)*T_matrix)*transpose(T_matrix)*transpose(x_2);       %Solving for coefficients matrix
p_coeff = fliplr(polyfit(t',x_noise2,PolDegree));
p = T_matrix*p_coeff';

plot(t,p,'Color',[1,.3,.5]);
legend('Displacement','w/Noise','Global Fit');
