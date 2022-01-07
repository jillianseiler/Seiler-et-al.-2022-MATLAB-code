function [controlFit] = controlFit (control, signal)

reg = polyfit(control, signal, 1); %least-squares linear regression fit control to signal

a = reg(1);
b = reg(2);

controlFit = a.*control + b; %save out new fitted control vector

