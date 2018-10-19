%PIECEWISE FUNCTION%
clear variables
close all
clc

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
        x(i) = abs(cos(pi/2*(t-tau1))).*exp(-(C*sin(2*pi*(t-tau1)))).*tan(pi./(1+exp(-(t-tau1))));
    elseif t > 0;
        x(i) = abs(cos(pi/2*(t+tau2))).*exp(-(C*sin(2*pi*(t+tau2)))).*tan(pi./(1+exp(-(t+tau2))));
    else
        x(i) = 0;
    end
    i = i+1;
end
t = linspace(-4,4,n);
x_noise = add_awgn_noise(x,snr); %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
plot(t,x,'Color',[.1,.1,.34]);
hold on;
plot(t,x_noise,'Color',[.2,.6,.15]);
hold all

%%
%PIECEWISE FIT%

t0 = -4;        %Beginning of range
tend = 4;       %End of range

n = 100;        %Total number of points within range
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
            x_3 = abs(cos(pi/2*(t(i)-tau1))).*exp(-(C*sin(2*pi*(t(i)-tau1)))).*tan(pi./(1+exp(-(t(i)-tau1))));     %Equation when t is negative
            x_noise3(i,1) = add_awgn_noise(x_3,snr);      %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
        elseif t(i) > 0;
            x_3 = abs(cos(pi/2*(t(i)+tau2))).*exp(-(C*sin(2*pi*(t(i)+tau2)))).*tan(pi./(1+exp(-(t(i)+tau2))));     %Equation when t is positive
            x_noise3(i,1) = add_awgn_noise(x_3,snr);      %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
        else
            x_3 = 0;       %Equation when t is zero
            x_noise3(i,1) = add_awgn_noise(x_3,snr);      %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
        end
    end
    
    T_matrix = [];
    for k = 1:length(t)-1
        T_matrix = [T_matrix t.^(k-1)];
    end
    
    p_coeff = fliplr(polyfit(t,x_noise3,k-1));
    p = T_matrix*p_coeff';
    
    %     a = inv(T_matrix'*T_matrix)*T_matrix'*x_3;
    %     poly = T_matrix*a;
    %     poly_final = [poly_final;poly];
    %     x_final = [x_final; x_3];
    t_final = [t_final t'];
    p_final = [p_final p'];
    
end

plot(t_final,p_final,'Color',[1,.3,.5]);
 legend('Displacement', 'w/Noise','Piecewise Fit');
