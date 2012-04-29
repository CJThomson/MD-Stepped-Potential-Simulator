clc
clear all
close all

x = [0.001:0.0001:3.0];
y = 4 * (x.^-12 - x.^-6);
plot(x,y,'k')
xlim([0 3.0])
ylim([-1.5 1.5])
xlabel('r/$\sigma$','interpreter','none')
ylabel('U(r)/$\epsilon$','interpreter','none')
print('ljPlot.tex','-depslatex');
