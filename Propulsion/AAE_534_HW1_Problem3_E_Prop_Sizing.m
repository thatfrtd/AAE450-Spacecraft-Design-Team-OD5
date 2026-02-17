%AAE_534_HW1_Problem_3
%% ____________________
%% INITIALIZATION

%Change These Values for what our system will use

m_u = 500;          %kg
T = 0.500;          %N
dv = 8000;          %m/s
a = 0.0001;         %kg/W
n = 0.5;
I_sp = linspace(25000, 32000, 1000);
g_0 = 9.80665;      %m/s^2
%% ____________________
%% CALCULATIONS
m_ppm_pr = ((a*T*g_0*I_sp)/(2*n)).*exp(dv./(I_sp*g_0))+m_u*(exp(dv./(I_sp*g_0))-1);
[min_m, idx] = min(m_ppm_pr);

min_ppm = (a*T*g_0*I_sp(idx))/(2*n)*exp(dv/(I_sp(idx)*g_0));
min_pr = m_u*(exp(dv/(I_sp(idx)*g_0))-1);

%% ____________________
%% OUTPUTS
figure(1)
plot(I_sp, m_ppm_pr, 'b', 'LineWidth', 1.5)
hold on
plot(I_sp(idx), min_m, 'r*')
hold off
xlabel('I_s_p^* (s)')
ylabel('m_p_p_u + m_p_r')
title('m_p_p_u + m_p_r vs I_s_p')
legend('', sprintf('min m_p_p_u + m_p_r, I_s_p = %.3f', I_sp(idx)))
grid on




