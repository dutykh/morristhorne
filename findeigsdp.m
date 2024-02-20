% The double precision version

clear
close all
format longG

M0 = readmatrix('data/M0_200.mat', 'FileType', 'text');
M1 = readmatrix('data/M1_200.mat', 'FileType', 'text');
M2 = readmatrix('data/M2_200.mat', 'FileType', 'text');

e = polyeig(M0, 1i*M1, M2);

figure(1);
plot(real(e), imag(e), 'bx'), grid off, hold off
xlabel('$\mathrm{Re}\Omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
ylabel('$\mathrm{Im}\Omega$', 'interpreter', 'LaTeX', 'fontsize', 14);
axis([-6 6 -5 1]);
title('Worm hole QNMs');
set(gcf, 'color', 'w');