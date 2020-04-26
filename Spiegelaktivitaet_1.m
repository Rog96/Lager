%% combine emg- and force data
% structure emg and force data
%   ch 1        ch 2        ...         ch n
%   Value 1     Value 1     Value 1     Value 1
%   Value 2     Value 2     Value 2     Value 2
%   ...         ...         ...         ...
%   Value n     Value n     Value n     Value n
dataSet = [emgData, forceData];
% number of channels
nCh = size(dataSet,2);

%% label
emgLabel={'L M. pect. major', 'R M. pect. major',...
    'L M. trapezius (Pars transversal)', 'R M. trapezius (Pars transversal)',...
    'L M. latissimus dorsi', 'R M. latissimus dorsi',...
    'L M. deltoideus (Pars acromialis)', 'R M. deltoideus (Pars acromialis)',...
    'L M. triceps brachii (Caput laterale)', 'R M. triceps brachii (Caput laterale)',...
    'L M. biceps brachii (Caput breve)', 'R M. biceps brachii (Caput breve)',...
    'L M. Flexor Carpi Radialis', 'R M. Flexor Carpi Radialis',...
    'L M. brachioradialis', 'R M. brachioradialis'};
forceLabel={strcat(char(hex2dec('03A3')),' Fz') 'R FD', 'L FD', char(hex2dec('03C6'))};

%% selected color
color1=[0.3010, 0.7450, 0.9330]; % blue
color2=[0.9290, 0.6940, 0.1250]; % orange

%% rectify emg signals
emgAbs = abs(emgData);

%% RMS EMG-Kanäle (windows 50 ms)
emgRms = sqrt(movmean(emgData.^2,floor(.05 *fS)));

%% fft 
dataSetFFT = zeros(length(dataSet),nCh);
for i=1:nCh
    dataSetFFT(:,i) = fft(dataSet(:,i));
    dataSetFFT(:,i) = dataSetFFT(:,i) / length(dataSet);
end
dataSetFFT = [dataSetFFT(1,:); 2*dataSetFFT(2:floor((length(dataSet))/2),:)];

% frequency axis
fAxis = transpose(linspace(0, fS/2, floor(length(dataSet)/2)));
fAxis = fAxis';

%% split in On-/Off phases
% determine number of phases (on / off)
nPhases = floor((max(M)-1)/2);
% parameter dataSet, emgAbs, emgRms
nParameters = 3;
% structure activity on respectively activity off
%   phase 1     phase 2     ...         phase n
%   dataset     dataset     ...         dataset
%   emgAbs      emgAbs      ...         emgAbs
%   emgRms      emgRms      ...         emgRms

activityOn = {nParameters,nPhases};
activityOff = {nParameters,nPhases};
nOn = zeros(1,nPhases);
nOff = zeros(1,nPhases);
t1Padding=fS;

for i=1:nPhases
    t1 = find(M > i*2-2,1);
    t2 = find(M > i*2-1,1);
    t3 = find(M > i*2,1);
    
    activityOn{1,i} = dataSet(t1-t1Padding:t2,:);
    activityOn{2,i} = emgAbs(t1-t1Padding:t2,:);
    activityOn{3,i} = emgRms(t1-t1Padding:t2,:);
    
    activityOff{1,i} = dataSet(t2:t3,:);
    activityOff{2,i} = emgAbs(t2:t3,:);
    activityOff{3,i} = emgRms(t2:t3,:);
    
    nOn(i) = length(activityOn{1,i});
    nOff(i) = length(activityOff{1,i});
end
% n for resampling
nOn = floor(mean(nOn));
nOff = floor(mean(nOff));
clearvars t1 t2 t3 i

%% resampling phases
tempOn = {nParameters,nPhases};
tempOff = {nParameters,nPhases};

for i=1:nParameters
    for j=1:nPhases
        tempOn{i,j} = resample(activityOn{i,j},nOn,length(activityOn{i,j}));
        tempOff{i,j} = resample(activityOff{i,j},nOff,length(activityOff{i,j}));
    end
end
% create new data structure
temp1On = zeros(nOn,nPhases);
temp1Off = zeros(nOff,nPhases);
activityOnNorm={nParameters,nCh};
activityOffNorm={nParameters,nCh};

% structure activityOnNorm respectively activityOffNorm
%   ch 1        ch 2        ...         ch n
%   dataset     dataset     ...         dataset
%   emgAbs      emgAbs      ...         emgAbs
%   emgRms      emgRms      ...         emgRms

for i=1:nParameters
    for j=1:size(tempOn{i,1},2)
        for k=1:nPhases
            temp1On(:,k) = tempOn{i,k}(:,j);
            temp1Off(:,k) = tempOff{i, k}(:,j);
        end
        activityOnNorm{i,j} =temp1On;
        activityOffNorm{i,j} =temp1Off;
    end
end
clearvars activityOn activityOff tempOn tempOff temp1On temp1Off i j k
%% mean curve from emgRms & force channels
dataSetNormOn=[activityOnNorm(3,1:size(emgData,2)),activityOnNorm(1,size(emgData,2)+1:nCh)];
dataSetNormOff=[activityOffNorm(3,1:size(emgData,2)),activityOffNorm(1,size(emgData,2)+1:nCh)];

outputOnCurve = zeros(nOn,nCh);
for i=1:nCh
    outputOnCurve(:,i) = mean(dataSetNormOn{1, i},2);
end

outputOffCurve = zeros(nOff,nCh);
for i=1:nCh
    outputOffCurve(:,i) = mean(dataSetNormOff{1, i},2);
end
outputCurve = [outputOnCurve; outputOffCurve];
toutputCurve = transpose(-t1Padding*1/fS:1/fS:(length(outputCurve)-t1Padding-1)*1/fS);
%%
chReference=17;
chSelect=7;

figure('Name','Measured Data',...
    'NumberTitle','off',...
    'Color','w',...
    'WindowState','maximized',...
    'units','normalized',...
    'outerposition',[0 0 1 1]);

annotation('rectangle',...
    'Color',[.8 .8 .8],...
    'units','normalize',...
    'position',[.01 .01 .98 .98]);

for i=1:6
subplot(3,2,i);
yyaxis left
plot(toutputCurve,outputCurve(:,i),...
    'color',color1,...
    'lineWidth',1);
emgYLim = roundn(max(outputCurve(:,i)), 1);
ylim([-10 emgYLim]);
yticks([0 .5*emgYLim emgYLim]);
ytickformat('%g µV');
yLeftLabel=ylabel(['{\color[rgb]{0.3010, 0.7450, 0.9330}',char(hex2dec('25A0')),...
    '\color[rgb]{0 0 0}','   \rightarrow RMS EMG (50 ms)}'],...
    'interpreter','tex');

yyaxis right
referenceForce = outputCurve(:,chReference) - mean(outputOffCurve(:,chReference));
plot(toutputCurve,referenceForce,...
    'color',color2,...
    'lineWidth',1.5);
xline(0,'--','color',[0,0,0]);
grid on;
forceYLim = roundn(max(referenceForce), 2);
ylim([-10 forceYLim]);
yticks([0 .5*forceYLim forceYLim]);
ytickformat('%g N');
yRightLabel=ylabel(['{\color[rgb]{0.9290, 0.6940, 0.1250}',char(hex2dec('25A0')),...
    '\color[rgb]{0 0 0}','   \rightarrow vertikale Bodenreaktionskraft}'],...
    'interpreter','tex');

xlim([min(toutputCurve),max(toutputCurve)]);
xtickformat('%g s');
xLabel = xlabel('\rightarrow Zeit');

hold on;
% legend with icons
markerSizeLegend = 7;
chlegend=plot(NaN,NaN,...
    'square',...
    'DisplayName',string(emgLabel(i)),...
    'MarkerSize',markerSizeLegend,...
    'MarkerEdgeColor',color1,...
    'MarkerFaceColor',color1);
refLegend=plot(NaN,NaN,...
    'square',...
    'DisplayName','Referenz',...
    'MarkerSize',markerSizeLegend,...
    'MarkerEdgeColor',color2,...
    'MarkerFaceColor',color2);

aL=legend([chlegend refLegend],'FontSize',9,...
    'Orientation','horizontal');
legend('boxoff')
aL.ItemTokenSize = [10;18];
pos_aL=get(aL, 'Position');
pos_aP = get(gca,'Position');
aL.Position = [pos_aP(1)+.5*pos_aL(3), pos_aP(2)+pos_aP(4)+.02, 0, 0];

% set position axis label
set(xLabel,'units','normalize')
xLabelPos=xLabel.Position;
xLabel.Position = [1,xLabelPos(2),0];
xLabel.HorizontalAlignment = 'right';

set(yLeftLabel,'units','normalize')
y1LabelPos=yLeftLabel.Position;
yLeftLabel.Position = [y1LabelPos(1),1,0];
yLeftLabel.HorizontalAlignment = 'right';

set(yRightLabel,'units','normalize')
y2LabelPos=yRightLabel.Position;
yRightLabel.Position = [y2LabelPos(1),1,0];
yRightLabel.HorizontalAlignment = 'right';

% Current axes or chart
ax = gca;
% set axis properties
ax.YAxis(1).Color = [.2,.2,.2];
ax.YAxis(2).Color = [.2,.2,.2];
ax.FontSize = 8;
ax.FontName = 'Helvetica';
end

%% mean/std activity on/off phases on the basis of rms
% structure outputParameter
%   ch 1    ch 2	... 	ch n
%   meanOn	meanOn	...     meanOn
%   stdOn	stdOn	...     stdOn
%   meanOff	meanOff	...     meanOff
%   stdOff	stdOff	...     stdOff
%   snr     snr     ...     snr

outputParameter = zeros(5,nCh);
phasePadding = fS;

for i=1:nCh
    outputParameter(1,i) = mean(dataSetNormOn{1, i}(phasePadding:end-phasePadding,:),'all');
    outputParameter(2,i) = std(dataSetNormOn{1, i}(phasePadding:end-phasePadding,:),0,'all');
    outputParameter(3,i) = mean(dataSetNormOff{1, i}(phasePadding:end-phasePadding,:),'all');
    outputParameter(4,i) = std(dataSetNormOff{1, i}(phasePadding:end-phasePadding,:),0,'all');
end
outputParameter(5,:) = outputParameter(1,:) ./ outputParameter(4,:);
clearvars dataSetParameterOn dataSetParameterOff i

%% fft activity ON phase
nOnPadding = nOn-2*phasePadding+1;
temp = zeros(nOnPadding,nPhases);
onFFT = zeros(floor(nOnPadding/2),nCh);

for i=1:nCh
    for j=1:nPhases
        temp(:,j) = fft(activityOnNorm{1, i}(phasePadding:size(activityOnNorm{1, i}) - phasePadding,j));
        temp(:,j) = temp(:,j) / nOnPadding;
    end
    temp1=mean(abs(temp),2);
    onFFT(:,i) = [temp1(1,:); 2*temp1(2:floor(nOnPadding/2),:)];
end

% frequency axis
fAxisOn = transpose(linspace(0, fS/2, floor(nOnPadding/2)));
clearvars temp temp1 i j

%% fft activity OFF phase
nOffPadding = nOff-2*phasePadding+1;
temp = zeros(nOffPadding,nPhases);
offFFT = zeros(floor(nOffPadding/2),nCh);

for i=1:nCh
    for j=1:nPhases
        temp(:,j) = fft(activityOffNorm{1, i}(phasePadding:size(activityOffNorm{1, i}) - phasePadding,j));
        temp(:,j) = temp(:,j) / nOffPadding;
    end
    temp1=mean(abs(temp),2);
    offFFT(:,i) = [temp1(1,:); 2*temp1(2:floor(nOffPadding/2),:)];
end

% frequency axis
fAxisOff = transpose(linspace(0, fS/2, floor(nOffPadding/2)));
clearvars temp temp1 i j

%%
temp=activityOnNorm{1, chSelect}(phasePadding:size(activityOnNorm{1, chSelect}) - phasePadding,:);
subplot(2,2,1)

pwelch(temp,[],[],[],fS)

[pxx,fPxx] = pwelch(temp,[],[],[],fS);
pxxMean = mean(pxx,2);
pxxMeanDb = pow2db(pxxMean);

energy = cumtrapz(fPxx,pxxMeanDb);
energyMax = find(energy == max(energy(:)));
if energyMax > 1
    fPxxEnd = fPxx(energyMax);
else
    fPxxEnd = 500;
end

[fPxxMean,~]=meanfreq(pxxMean,fPxx,[0 fPxxEnd]);
% fPxxMean = sum(pxxMean(fPxx<fPxxEnd,:).*fPxx(fPxx<fPxxEnd,:))/sum(pxxMean(fPxx<fPxxEnd,:))

subplot(2,2,2)
plot(fPxx,pow2db(pxx));hold on
yline(0,'k');

xlim([0 fS/2])
grid on
subplot(2,2,3)
plot(fPxx,pxxMean)
xline(fPxxMean,'-.','color',[0,0.4470, 0.7410]);
xlim([0 fS/2])
grid on

subplot(2,2,4)
plot(fPxx,pow2db(pxxMean))
ylim([-30 50])
yline(0,'k');
xline(fPxxMean,'--','color',[0,0.4470, 0.7410]);
rectangle('Position',[fPxxEnd -100 fS/2 200],'FaceColor',[.7 .7 .7 .25],...
   'EdgeColor',[1 1 1])
xlim([0 fS/2])
grid on
hold on
area(fPxx(1:energyMax),pow2db(pxxMean(1:energyMax)),...
    'FaceColor',[0, 0.4470, 0.7410],'EdgeColor',[0, 0.4470, 0.7410],'FaceAlpha',0.2)

set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',8);
hold off;
%%
s=spectrogram(emgData(:,10));
spectrogram(emgData(:,10),128,28,128,fS,'yaxis','MinThreshold',-40)
colormap Summer

%% temp
ecgChannels=[emgData(:,1), emgData(:,2), emgData(:,5), emgData(:,6)];
ecgThreshold=2.5*std(abs(ecgChannels(:,1)));

%% ECG
% find R-peak from QRS-complex
[ecgValue,ecgPos] = findpeaks(ecgChannels(:,1),'MinPeakDistance',0.4*fS,...
    'MinPeakHeight',ecgThreshold);
% calculate Heartrate from EKG
heartRate = 60 ./ (1/fS * diff(ecgPos));
% find outlier and correct
[correctHeartRate,posErrorHeartRate] = rmoutliers(heartRate,'median');
nECGErrorMarker = find(posErrorHeartRate == 0);
correctPosEKG =ecgPos(nECGErrorMarker,:);

plot(t,ecgChannels(:,1),t(ecgPos),ecgValue,'or');hold on;
plot(t([1 end]),[1 1] * ecgThreshold,'k--');hold on;


%% ECG template
%number of arguments
nECGTemplate = floor(fS * 60 / mean(correctHeartRate));

%mean template from selected channel
%temp0=zeros(nECGTemplate,length(correctPosEKG) - 1);
temp0=zeros(nECGTemplate,5);
% for i = 1:length(correctPosEKG)-1
sonne=dataSet(:,1);
for i = 1:5
    temp1 = sonne(ecgPos(i):ecgPos(i+1),1);
    temp0(:,i) = resample(temp1,nECGTemplate, length(temp1));
end
ecgTemplate = mean(temp0,2);
%clearvars temp0 temp1

%%
emgTemp = emgData(1:length(visualECG),1);
emgTempFFT = fft (emgTemp); % Fouriertransformation
emgTempFFT = emgTempFFT / length(emgTemp);
emgTempFFT = [emgTempFFT(1,:); 2*emgTempFFT(2:floor((length(emgTemp))/2),:)];
fAxisEmgTempFFT = transpose(linspace(0, fS/2, floor(length(emgTemp)/2)));
fAxisEmgTempFFT = fAxisEmgTempFFT';

emgTempT = transpose(0:1/fS:(length(emgTemp)-1)*1/fS);
% Plot
subplot(3,2,1)
plot(emgTempT,emgTemp);
subplot(3,2,2)
plot(fAxisEmgTempFFT,abs(emgTempFFT));

%%
ecgTemp = ecgTemplate;
ecgTempFFT = fft (ecgTemp); % Fouriertransformation
ecgTempFFT = ecgTempFFT / length(ecgTemp);
ecgTempFFT = [ecgTempFFT(1,:); 2*ecgTempFFT(2:floor((length(ecgTemp))/2),:)];

ecgTempFFT = resample(ecgTempFFT,length(emgTempFFT),length(ecgTempFFT));

fAxisEcgTempFFT = transpose(linspace(0, fS/2, floor(length(emgTemp)/2)));
fAxisEcgTempFFT = fAxisEcgTempFFT';

ecgTempT = transpose(0:1/fS:(length(ecgTemp)-1)*1/fS);
% Plot
subplot(3,2,3)
plot(ecgTempT,ecgTemp);
subplot(3,2,4)
plot(fAxisEcgTempFFT,abs(ecgTempFFT));


%%
bla1=emgTempFFT;
bla2=ecgTempFFT;

bla3=abs(bla2)-abs(bla1);
bla4=angle(bla2)-angle(bla1);

bla5 = bla3.*cos(bla4);
bla6= bla3 .* sin(bla4);
bla7= complex(bla5,bla6);
bla8=bla1-bla2;

% 
Output=ifft(bla7);

OutputFFT = fft (Output); % Fouriertransformation
OutputFFT = OutputFFT / length(Output);
OutputFFT = [OutputFFT(1,:); 2*OutputFFT(2:floor((length(Output))/2),:)];

fAxisOutputFFT = transpose(linspace(0, fS/2, floor(length(Output)/2)));
fAxisOutputFFT = fAxisOutputFFT';

OutputT = transpose(0:1/fS:(length(Output)-1)*1/fS);
% Plot
subplot(3,2,5)
plot(OutputT,Output);
subplot(3,2,6)
plot(fAxisOutputFFT,abs(OutputFFT));




% visualECG = repmat(ecgTemplate,5,1);
% tvisualECG = 0:1/fS:(length(visualECG)-1)*1/fS;
% tvisualECG = tvisualECG';

%%
[c,l] = wavedec(emgTemp,3,'db2');
approx = appcoef(c,l,'db2');
[cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
wav1 = waverec(c,l,'db2');
subplot(4,1,1)
plot(approx)
title('Approximation Coefficients')
subplot(4,1,2)
plot(cd3)
title('Level 3 Detail Coefficients')
subplot(4,1,3)
plot(cd2)
title('Level 2 Detail Coefficients')
subplot(4,1,4)
plot(cd1)
title('Level 1 Detail Coefficients')


%%
wname = 'db8';
lev = 5;
[c,l] = wavedec(emgTemp,lev,wname);
det = wrcoef('d',c,l,wname,2);
figure
plot(emgTemp);hold on;
plot(det)
title('Level 1 Detail')
%%
level = 6;
wpt = wpdec(emgTemp,level,'sym8');
[Spec,Time,Freq] = wpspectrum(wpt,fS,'plot');
figure;
windowsize = 128;
window = hanning(windowsize);
nfft = windowsize;
noverlap = windowsize-1;
[S,F,T] = spectrogram(emgTemp,window,noverlap,nfft,fS);
imagesc(T,F,log10(abs(S)))
set(gca,'YDir','Normal')
xlabel('Time (secs)')
ylabel('Freq (Hz)')
title('Short-time Fourier Transform spectrum')
%%
[wt,f] = cwt(emgTemp,1500);
xrec = icwt(wt,f,[80 500],'SignalMean',mean(emgTemp));
subplot(2,1,1)
plot(emgTemp)
grid on
title('Original Data')
subplot(2,1,2)
plot(xrec)
grid on
title('Bandpass Filtered Reconstruction [0.030 0.070] Hz')
%%
[cA1,cD1] = dwt(emgTemp,'db1');
l_s = length(emgTemp);
A1 = upcoef('a',cA1,'db3',1,l_s);
D1 = upcoef('d',cD1,'db3',1,l_s);
subplot(1,2,1); plot(A1); title('Approximation A1')
subplot(1,2,2); plot(D1); title('Detail D1')
A0 = idwt(cA1,cD1,'db1',l_s);
err = max(abs(emgTemp-A0));
[C,L] = wavedec(emgTemp,3,'db1');
cA3 = appcoef(C,L,'db1',3);
cD3 = detcoef(C,L,3);
cD2 = detcoef(C,L,2);
cD1 = detcoef(C,L,1);
A3 = wrcoef('a',C,L,'db1',3);
D1 = wrcoef('d',C,L,'db1',1); 
D2 = wrcoef('d',C,L,'db1',2); 
D3 = wrcoef('d',C,L,'db1',3);
subplot(2,2,1); plot(A3);  
title('Approximation A3') 
subplot(2,2,2); plot(D1);  
title('Detail D1') 
subplot(2,2,3); plot(D2);  
title('Detail D2') 
subplot(2,2,4); plot(D3);  
title('Detail D3')
A0 = waverec(C,L,'db1');  
err = max(abs(emgTemp-A0)) 
% subplot(2,1,1);plot(emgTemp);title('Original'); axis off 
% subplot(2,1,2);plot(A3);title('Level 3 Approximation'); 
% axis off
%%
emgTemp = emgData(:,5);
 wtecg = modwt(emgTemp,8,'db4');  % obtain the MODWT at three levels of resolution -- scales, 2,4,8
  dt = 1/1500; % data sampled at 180 Hz
 t = 0:dt:(length(emgTemp)*dt)-dt;
 ecgmra = modwtmra(wtecg,'haar');
 subplot(10,1,1);
 plot(t,emgTemp); title('Original Data');
 for kk = 1:size(ecgmra,1)
   subplot(10,1,kk+1)
   plot(t,ecgmra(kk,:));
 end
figure;
 ts = sum(ecgmra,1);
 plot(t,[emgTemp ts'])
 grid on;
 %%
%  qq1= [emgTemp,ecgmra'];
ecgmra1=ecgmra(1:9,:);
 ecgChannels1 = zscore(ecgmra1,0,'all');
 %ecgChannels1=ecgChannels1';
 [Zica, W, T, mu] = fastICA(ecgChannels1,2);
 qqBla=std(emgTemp);
 qqBla1=qqBla*Zica;
 qqBla2=qqBla1';
 %%
 zzz=Zica';

%% Filter emg channels (highpass) [remove ecg signal]
% Stopband Frequency, Passband Frequency, Stopband Attenuation (dB), Passband Ripple (dB)
filterPara =[80, 85, 80, 1];
h = fdesign.highpass('fst,fp,ast,ap', filterPara(1), filterPara(2),...
    filterPara(3), filterPara(4), fS);
Hd = design(h, 'cheby2', 'MatchExactly', 'stopband', 'SOSScaleNorm', 'Linf');

emgFil = zeros(length(dataSet),nCh);
for i=1:nCh
    emgFil(:,i) = filter(Hd,dataSet(:,i));
end
clearvars filterPara h Hd








