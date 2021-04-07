maxMotorTorque = 80; % in Nm
maxRegenTorque = 80; % in Nm
maxMotorTemp = 135;   % in Fahrenheit
motorVoltage = 140;
motorTorques = [80 80 80 80 80 80 80 80 71.62 63.66 57.30 52.09 47.75 44.07 40.93];
motorSpeeds = [0 500 1000 1500 2000 2500 3000 3580.977845 4000 4500 5000 5500 6000 6500 7000] * pi/30;

plot(motorSpeeds, motorTorques);
axis([0, 1000, 0, 100]);

maxEngineTorque = 80; % in Nm

maxEngineTemp = 240;  % in Fahrenheit


engineTorqueMap = [0,0;4,4.1;8,8.1;10.7,10.8;9.2,9.3;11.8,11.9;14.8,14.9;17.3,17.4;
    18,18.1;19.2,19.3;18,18.1;17.8,17.9;17.5,17.6;16,16.1;15,15.1] * 1.3558;
engineTorqueBreakpoints = [0:1.378:19.3];
engineSpeedBreakpoints = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14] * 1000;


initialFuel = 100; % NEED TO SET TO CORRECT AMOUNT
