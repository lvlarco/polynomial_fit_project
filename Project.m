%PROJECT%
clear variables
close all
clc

%VARIABLES%
C = 0.9;
tau1 = 0.4;
tau2 = 0.15;

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
plot(t,x);
hold on

%%
%GLOBAL FIT%

PolDegree = 66;         %Degree of polynomial
n = PolDegree + 1;      %Number of points in fit
i = 1;

%For-loop discretizes equation and stores values in vector x
for t = linspace(-4,4,n);       %For-loop control statement
    if t<0;
        x_2(i) = abs(cos(pi/2*(t-tau1))).*exp(-(C*sin(2*pi*(t-tau1)))).*tan(pi./(1+exp(-(t-tau1))));        %Equation when t is negative
    elseif t > 0;
        x_2(i) = abs(cos(pi/2*(t+tau2))).*exp(-(C*sin(2*pi*(t+tau2)))).*tan(pi./(1+exp(-(t+tau2))));        %Equation when t is positive
    else 
        x_2(i) = 0;     %Equation when t is zero
    end
    i = i+1;
end

%Approach to solve polynomial fit: a = (T^t*T^-1)^(-1)*T^t*x
t = linspace(-4,4,n);
T_matrix = fliplr(vander(t));       %Creting T matrix
a = inv(transpose(T_matrix)*T_matrix)*transpose(T_matrix)*transpose(x_2);       %Solving for coefficients matrix
poly = T_matrix*a;      %Polynomial vector

plot(t,poly,'Color',[1,.3,.5]);
ylabel('X');xlabel('t');
legend('Displacement','Global Fit');
