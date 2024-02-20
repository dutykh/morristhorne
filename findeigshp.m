% The multi-precision version

clear
close all
format longG

addpath('../advanpix/'); % please set the path to advanpix repository

mp.Digits(200);

M0 = mp.read('M0_200.mat');
M1 = mp.read('M1_200.mat');
M2 = mp.read('M2_200.mat');

e = polyeig(M0, mp('1i')*M1, M2);

figure(1);
plot(real(e), imag(e), 'bx'), grid off, hold off
xlabel('$\mathrm{Re}\Omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
ylabel('$\mathrm{Im}\Omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
axis([-6 6 -5 1]);
title('Worm hole QNMs');
set(gcf, 'color', 'w');