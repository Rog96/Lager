clearvars; close all; clc; 
%% Import files (*.m)
% directory & trial
% exportIndex is trial*_EMG and trial*_Force
Directory = 'C:\Users\klaus\Desktop\Pilotexperiment_Spiegelaktivitaet\Proband_01\';
trial = 'P01_08';
% sample frequency
fS = 1500;
% calibration factor forceplates 
KF = 0.2863;

%% import EMG data
load(cat(2,Directory,trial, '_EMG.mat'),'Data');
% measured values ??number, number of channels
emgCh = [length(Data{1,1}),length(Data)];
% convert cellArray to array
emgData = zeros(emgCh);2 #hier wurde ein leeres Array mit Nullen kreiert was für die Datenauswertung genutzt wird 
for i=1:emgCh(2) #i=1 entspricht der Zahl die für die Länge der Zeilen zugeteilt wurde (was in '15' definiert wurde) emgCH 2 ist die Zahl zwei die für die Anzahl der Spalten steht  
    emgData(:,i) = Data{1,i}; 
end
clearvars emgCh Data;

%% import Force data
load(cat(2,Directory,trial, '_Force.mat'),'Data','Markers');
% measured values ??number, number of channels
forceCh = [length(Data{1,1}),length(Data)];
% convert cellArray to array
forceData = zeros(forceCh);
for i=1:forceCh(2)
    forceData(:,i) = Data{1,i};
end
clearvars forceCh Data i

%% sync files #beide Datensätze auf die gleiche Länge bringen und einen Startpunkt festlegen 
emgSync = find(emgData(:,end) > 0, 1);
forceSync = find(forceData(:,end) > 0, 1);

if emgSync < forceSync
    forceData = forceData(forceSync-emgSync+1:end,:);
else
    emgData = emgData(emgSync-forceSync+1:end,:);
end
forceData = forceData(1:min(length(emgData),length(forceData)),:);
emgData = emgData(1:min(length(emgData),length(forceData)),:);
clearvars emgSync forceSync

%% sync Markers
M = zeros(length(emgData),1);
for i =1:length(Markers)
    temp = find(forceData(:,1) >= Markers(i),1);
    M(temp,1)=i;
end
emgData = emgData(:,2:end -1);
forceData = forceData(:,2:end -1);
%% detrend emg channels
emgData = detrend(emgData);
clearvars Markers temp i

%% calibration & offset force channels
groundReactionForce = (forceData(:,1) + forceData(:,2)) * KF;
forceData(:,1) = groundReactionForce;
forceData(:,2) = [];
clearvars groundReactionForce offset

%% create dataset
t = transpose(0:1/fS:(length(emgData)-1)*1/fS);
dataSet = [];
dataSet(:,1) = t;
dataSet(:,2) = M;
dataSet = [dataSet,emgData,forceData];

%% export dataset and formating output
fileName = cat(2,Directory,trial,'.txt');
export = dataSet;
export(:,1) = round(export(:,1),4);
export(:,3:end) = round(export(:,3:end),1);

writematrix(export,fileName,'Delimiter','tab')
clearvars KF fileName trial Directory export dataSet

