%panelAttitude = simOut.yout{5}.Values;
%umbra = simOut.yout{8}.Values.Data;
%throttle = simOut.yout{9}.Values.Data;
%relAngle = simOut.yout{10}.Values.Data;
%battState = simOut.yout{11}.Values.Data;
%chargeEfficiency = simOut.yout{12}.Values.Data;
%powerGen = simOut.yout{13}.Values.Data;
%optAngle = simOut.yout{14}.Values.Data;
%panelRotationAngle = simOut.yout{6}.Values.Data;
%times

x = readmatrix("optimalBodyAngle.csv");

endTime = 365.25*24*3600;
t = linspace(0, endTime, length(x));

plot(t/3600/24, x(:,2))
xlabel('Time (days)');
ylabel('Body Rotation Angle (degrees)');
title('Optimal Body Rotational Angle Over a Year');
xlim([0,365])
