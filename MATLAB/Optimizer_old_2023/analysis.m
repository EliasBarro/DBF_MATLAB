%https://github.com/Georacer/ardupilog

% Open file prompt to select *.bin file, then convert to matlab struct
log = Ardupilog('rawData.BIN').getStruct;

%Flight path
plot3(log.GPS.Lng, log.GPS.Lat, log.GPS.Alt)

%Accel
time = log.IMU_1.TimeS;
plot(time, log.IMU_1.AccX, time, log.IMU_1.AccY, time, log.IMU_1.AccZ)