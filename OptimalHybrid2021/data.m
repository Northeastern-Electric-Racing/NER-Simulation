maxMotorTorque = 80; % in Nm
maxRegenTorque = 80; % in Nm
maxMotorTemp = 135;   % in Fahrenheit
motorVoltage = 140;
motorTorques = [80 80 80 80 80 80 80 80 71.62 63.66 57.30 52.09 47.75 44.07 40.93];
motorSpeeds = [0 500 1000 1500 2000 2500 3000 3580.977845 4000 4500 5000 5500 6000 6500 7000] * pi/30;

figure(1);
plot(motorSpeeds, motorTorques);
axis([0, 1000, 0, 100]);

maxEngineTorque = 26.16; % in Nm

maxEngineTemp = 240;  % in Fahrenheit

efficiencyTrq = [0:5.714:80];

efficiency = [
    52 78 80 82 82 82 82 82 82 80 78 78 76 74 70; 
    52 78 80 82 82 82 82 82 82 80 78 78 76 74 70;
    52 78 84 86 86 86 86 86 86 86 86 84 84 82 82;
    52 80 86 86 88 88 88 88 88 88 88 86 86 84 84;
    52 80 86 88 88 88 88 88 88 88 88 88 86 86 86;
    52 80 86 88 90 90 90 90 90 90 90 90 88 88 88;
    52 82 88 90 91 91 92 92 92 92 92 92 92 91 91;
    52 82 88 91 91 92 93 93 93 93 93 93 93 93 93;
    52 84 90 92 92 93 93 93 93 93 93 93 93 93 93; 
    62 86 91 93 93 93 93 93 93 93 93 93 93 93 93;
    66 86 91 93 93 93 93 93 93 93 93 93 93 93 93; 
    70 86 92 93 93 93 93 93 93 93 93 93 93 93 93; 
    70 86 92 93 93 93 93 93 93 93 93 93 93 93 93;
    72 86 91 93 93 93 93 93 93 93 93 93 93 93 93;
    72 86 91 93 93 93 93 93 93 93 93 93 93 93 93];


engineTorqueMap = [0,0;4,4.1;8,8.1;10.7,10.8;9.2,9.3;11.8,11.9;14.8,14.9;17.3,17.4;
    18,18.1;19.2,19.3;18,18.1;17.8,17.9;17.5,17.6;16,16.1;15,15.1] * 1.3558;
engineTorqueBreakpoints = [0:1.378:19.3];
engineSpeedBreakpoints = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14] * 1000;

figure(2);
plot(engineSpeedBreakpoints, engineTorqueBreakpoints);

initialFuel = 2.834; % in Liters 
