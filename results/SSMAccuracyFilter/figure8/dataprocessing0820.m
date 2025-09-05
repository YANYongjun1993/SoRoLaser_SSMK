%% load new control data
clearvars;
clc
clear all
close all
% Load new controller logsout with random inputs
load 202507102123.mat
sensor1_Quater1 = get(logsout,"sensor1_Quater1");
time = sensor1_Quater1.Values.Time;
Sensor1Quater1 = sensor1_Quater1.Values.data;
sensor1_Quater2 = get(logsout,"sensor1_Quater2");
Sensor1Quater2 = sensor1_Quater2.Values.data;
sensor1_Quater3 = get(logsout,"sensor1_Quater3");
Sensor1Quater3 = sensor1_Quater3.Values.data;
sensor1_Quater4 = get(logsout,"sensor1_Quater4");
Sensor1Quater4 = sensor1_Quater4.Values.data;
sensor1_x = get(logsout,"sensor1_x");
Sensor1X = sensor1_x.Values.data;
sensor1_y = get(logsout,"sensor1_y");
Sensor1Y = sensor1_y.Values.data;
sensor1_z = get(logsout,"sensor1_z");
Sensor1Z = sensor1_z.Values.data;
Sensor1 = [Sensor1X';Sensor1Y';Sensor1Z';Sensor1Quater1';Sensor1Quater2';Sensor1Quater3';Sensor1Quater4'];
sensor2_Quater1 = get(logsout,"sensor2_Quater1");
Sensor2Quater1 = sensor2_Quater1.Values.data;
sensor2_Quater2 = get(logsout,"sensor2_Quater2");
Sensor2Quater2 = sensor2_Quater2.Values.data;
sensor2_Quater3 = get(logsout,"sensor2_Quater3");
Sensor2Quater3 = sensor2_Quater3.Values.data;
sensor2_Quater4 = get(logsout,"sensor2_Quater4");
Sensor2Quater4 = sensor2_Quater4.Values.data;
sensor2_x = get(logsout,"sensor2_x");
Sensor2X = sensor2_x.Values.data;
sensor2_y = get(logsout,"sensor2_y");
Sensor2Y = sensor2_y.Values.data;
sensor2_z = get(logsout,"sensor2_z");
Sensor2Z = sensor2_z.Values.data;
% concatenating the position and quaterion
Sensor2 = [Sensor2X';Sensor2Y';Sensor2Z';Sensor2Quater1';Sensor2Quater2';Sensor2Quater3';Sensor2Quater4'];

flag = get(logsout,"flag");
flagStatus = flag.Values.data;

Motor1_pos = get(logsout,"Motor1_pos_req_qc_mode2");
Motor1Pos = Motor1_pos.Values.data;
Motor2_pos = get(logsout,"Motor2_pos_req_qc_mode2");
Motor2Pos = Motor2_pos.Values.data;
Motor3_pos = get(logsout,"Motor3_pos_req_qc_mode2");
Motor3Pos = Motor3_pos.Values.data;
Motor4_pos = get(logsout,"Motor4_pos_req_qc_mode2");
Motor4Pos = Motor4_pos.Values.data;
uMotorPos = [Motor1Pos';Motor2Pos';Motor3Pos';Motor4Pos'];

yMotor = get(logsout,"yMotor");
yMotorCmd = yMotor.Values.data;
zMotor = get(logsout,"zMotor");
zMotorCmd = zMotor.Values.data;
initialCmd = [yMotorCmd';zMotorCmd'];

Cable1_tension = get(logsout,"Cable_tension_1");
Cable1Tension = Cable1_tension.Values.data;
Cable2_tension = get(logsout,"Cable_tension_2");
Cable2Tension = Cable2_tension.Values.data;
Cable3_tension = get(logsout,"Cable_tension_3");
Cable3Tension = Cable3_tension.Values.data;
Cable4_tension = get(logsout,"Cable_tension_4");
Cable4Tension = Cable4_tension.Values.data;
cableTension = [Cable1Tension';Cable2Tension';Cable3Tension';Cable4Tension'];
%% setup matrix
zu = {};
% spdRegFactor = 1
% keptIdxs = 4312:9339;
% y_eq = Sensor1(:,9468);
% spdRegFactor = 5
% keptIdxs = 4629:5629;
% y_eq = Sensor1(:,5768);
% spdRegFactor = 10
% keptIdxs = 2037:2554;
% y_eq = Sensor1(:,2733);
% spdRegFactor = 15
% keptIdxs = 3254:3604;
% y_eq = Sensor1(:,3802);
% spdRegFactor = 20
keptIdxs = 2280:2517;
y_eq = Sensor1(:,2788);
zu{1,1} = time(keptIdxs)'; 
zu{1,2} = (Sensor1(:, keptIdxs) - y_eq);
zu{1,3} = uMotorPos(:, keptIdxs); 
zu{1,4} = initialCmd(:, keptIdxs); 
zu{1,5} = cableTension(:, keptIdxs); 

save("Figure8BasicSpdRegFactor20.mat","zu")
