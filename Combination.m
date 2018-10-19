%PROJECT WITHOUT NOISE%
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
plot(t,x,'red', 'LineWidth',3);
hold all

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

%Plotting global fit
plot(t,poly,'green');
hold all;

%%
%PIECEWISE FIT%

t0 = -4;        %Beginning of range
tend = 4;       %End of range

n = 175;        %Total number of points within range
window_size = 5;        %Number of points in the window
number_of_windows = ceil(n/(window_size-1));        %Number of windows in range

delta_t = tend - t0;
step_size = delta_t/n;      %Constant step size throughout range

poly_final = [];
t_final = [];
x_final = [];
p_final =[];

for j = 1:number_of_windows      %From 1 to number of windows within range
    
    start_t = (j-1)*step_size*(window_size-1)+t0;        %Beginning value for t at each window
    end_t = start_t+step_size*(window_size-1);       %End value of t at each window
    t = start_t :step_size:end_t;
    t = t';
    
    for i = 1:length(t);
        if t(i) < 0;
            x_3(i,1) = abs(cos(pi/2*(t(i)-tau1))).*exp(-(C*sin(2*pi*(t(i)-tau1)))).*tan(pi./(1+exp(-(t(i)-tau1))));
        elseif t(i) > 0;
            x_3(i,1) = abs(cos(pi/2*(t(i)+tau2))).*exp(-(C*sin(2*pi*(t(i)+tau2)))).*tan(pi./(1+exp(-(t(i)+tau2))));
        else
            x_3(i,1) = 0;
        end
    end
    
    T_matrix = [];
    for k = 1:length(t)-1       %Creating Vandermonde matrix
        T_matrix = [T_matrix t.^(k-1)];
    end
    
    p_coeff = fliplr(polyfit(t,x_3,k-1));       %Calculating polynomial coefficient using polyfit function
    p = T_matrix*p_coeff';      %Calculating polynomial fit vector
    
    t_final = [t_final t'];     %Storing t and p data
    p_final = [p_final p'];
    
end

%Plotting piecewise fit
plot(t_final,p_final,'magenta');
hold all;


%%
%OPTIMAL FIT%
start =  -4;
endd = 4;
start_pt = start;        %Beginning of range
last_pt = endd;       %End of range

n = 300;        %Total number of points within range

delta_t = last_pt - start;
step_size = delta_t/n;      %Constant step size throughout range

t_final = [];
p_final =[];
count = 0;
windowss = [];
while start_pt < last_pt
    window_size = 1;
    Rsq = 1;
    x_4 = [];
    p = [];
    while Rsq >= 0.9999     %Criteria to decide window and polynomial size. R-square dependent
        window_size = window_size + 2;      %Window must be odd
        end_pt = start_pt+step_size*(window_size-1);
        t = start_pt:step_size:end_pt;      %Setting up time vector
        t = t';
        
        for i = 1:length(t);        %Discretizing equation
            if t(i) < 0;
                x_4(i,1) = abs(cos(pi/2*(t(i)-tau1))).*exp(-(C*sin(2*pi*(t(i)-tau1)))).*tan(pi./(1+exp(-(t(i)-tau1))));
            elseif t(i) > 0;
                x_4(i,1) = abs(cos(pi/2*(t(i)+tau2))).*exp(-(C*sin(2*pi*(t(i)+tau2)))).*tan(pi./(1+exp(-(t(i)+tau2))));
            else
                x_4(i,1) = 0;
            end
        end
        T_matrix = [];      %For-loop creating Vandermond matrix
        for k = 1:length(t)-1
            T_matrix = [T_matrix t.^(k-1)];
        end
        
        p_coeff = fliplr(polyfit(t,x_4,k-1));       %Calculating coefficients for polynomial
        p = T_matrix*p_coeff';      %Polynomial fit vector
        
        p_resid = x_4 - p;       %Residual between original data and polynomial fit
        SS_resid = sum(p_resid .^2);        %Residual sum of squares
        SS_total = (length(x_4)-1) * var(x_4);      %Multiplying by variance
        Rsq = 1 - SS_resid/SS_total;        %Calculating R squared value
        Rsq_adj = Rsq * (length(x_4)-1)/(length(x_4)-length(p_final));      %Adjusted R squared value
        plot(t_final, p_final);
    end
    t_final = [t_final t'];     %Creating total matrix for t
    p_final = [p_final p'];     %Creating total polynomial fit matrix
    windowss = [windowss window_size];      %Storing the size of each window
    start_pt = end_pt;
    count = count + 1;      %Counts how many windows are being created
    
end

%Plotting optimal fit and showing window count
disp(count)
plot(t_final, p_final, 'b');
legend('Displacement','Global Fit','Piecewise Fit','Optimal Fit');
xlim([start endd])