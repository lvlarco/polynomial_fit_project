%OPTIMAL FUNCTION%
clear variables
close all
clc

%Defining parameters
C = 0.9;
tau1 = 0.4;
tau2 = 0.15;
snr = 15;

%%
%DISCRETIZING AND PLOTTING THE DISPLACEMENT EQUATION%
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
plot(t,x,'green','LineWidth',2);
hold on;
plot(t,x_noise,'red');
hold all

%%
%OPTIMAL FIT%
start =  -4;
endd = 4;
t0 = -4;
start_pt = start;        %Beginning of range
last_pt = endd;       %End of range

n = 300;        %Total number of points within range 275

delta_t = last_pt - t0;
step_size = delta_t/n;      %Constant step size throughout range
t_final = [];
p_final =[];
count = 0;
windowss = [];
while start_pt < last_pt        %While-loop defining range
    window_size = 1;
    Rsq = 1;
    x_4 = [];
    p = [];
    while Rsq >= 0.9999     %R-square criteria to define window size
        window_size = window_size + 2;
        end_pt = start_pt+step_size*(window_size-1);
        t = start_pt:step_size:end_pt;
        t = t';
        
        for i = 1:length(t);
            if t(i) < 0;
                x_4 = abs(cos(pi/2*(t(i)-tau1))).*exp(-(C*sin(2*pi*(t(i)-tau1)))).*tan(pi./(1+exp(-(t(i)-tau1))));
                x_noise4(i,1) = add_awgn_noise(x_4,snr);      %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
            elseif t(i) > 0;
                x_4 = abs(cos(pi/2*(t(i)+tau2))).*exp(-(C*sin(2*pi*(t(i)+tau2)))).*tan(pi./(1+exp(-(t(i)+tau2))));
                x_noise4(i,1) = add_awgn_noise(x_4,snr);      %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
            else
                x_4 = 0;
                x_noise4(i,1) = add_awgn_noise(x_4,snr);      %Adding noise - add_awgn_noise function created by Mathuranathan Viswanathan
            end
        end
        T_matrix = [];      %For-loop creating Vandermond matrix
        for k = 1:length(t)-1
            T_matrix = [T_matrix t.^(k-1)];
        end
        
        p_coeff = fliplr(polyfit(t,x_noise4,k-1));       %Calculating coefficients for polynomial
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
    start_pt = end_pt+step_size;
    count = count + 1;      %Counts how many windows are being created
    
end
%Displays total count and plots optimal fit
disp(count)
plot(t_final, p_final, 'blue');
ylabel('X');xlabel('t');
legend('Displacement','w/Noise','Optimal Fit');

xlim([start endd])


