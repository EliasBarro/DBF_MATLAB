function crit_info(x_hist)

x=x_hist(1,3);
increment = 1;
    
while(x_hist(increment,3)==x)
    increment = increment + 1;
end

print=["takeoff distance (m): ", x_hist(increment,14)];
disp(print);
print3 = ["takeoff distance (ft)" ,x_hist(increment,14)*3.28084];
disp(print3)

print2=["max cruise speed: ",max(x_hist(:,2))];
disp(print2);