% intermediate functions

%% center clipping

k = -32:0.2:32;
v12 = sin(2*pi*k/32);
y12 = cclipper(v12,0.1);
subplot(2,1,1)
plot(k, v12)
ylabel('v12[k]')
grid on
subplot(2,1,2)
plot(k, y12)
ylabel('y12[k]')
grid on
%%
