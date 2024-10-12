%% Bugfixing tutorial


%% plotting
x = 1:.01:10;
y = x.^3 - x^2 + 2;
str = 5;

z = y/x;

figure
subplot(1,2,3);
plot(x,y, 'x');
xlim([10 -5]);
ylim([-10 1i]);
title(str);
subplot(2,2,3);
plot(z, x);




%% matrix math

A = 3*ones(10) - eye(5); 
B = 1:10;
C = A*B;

D = A\B;

for i = 1:size(B)
    num(i) = str2double(B(ii));
    
    
end


disp(num);



    






