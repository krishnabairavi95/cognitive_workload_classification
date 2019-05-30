%% Load data
clc
clear all
close all
mat = dir('*.mat');  
for q = 1:length(mat) 
    eeg(q)=load(mat(q).name); 
end
%% Setting variables
%4.4.2: To extract the EEG signal for the frequency range (0.01-40Hz) and apply filter
fs=512; % sampling frequency
T=1/fs; % time period
%Delta waves (0.5Hz to 3 Hz)
%Theta waves (3Hz to 8 Hz)
%Alpha waves (8Hz to 13 Hz)
%Beta waves (13 Hz to 38Hz)
%Gamma waves (38 Hz to 42 Hz)
freq=0:(fs/length(eeg(1).Data.EEG.raw.value)):40;

%% Filter for EEG data frequency range
w1=0.01/(0.5*fs);   %lower cut off frequency for EEG waves
w2=42/(0.5*fs);     %upper cut off frequency for EEG waves
wn=[w1,w2];         %setting filter coefficients
[b,a]=butter(2,wn); %creating filter
for i=1:length(mat)
    feeg(i,:)=filter(b,a,eeg(i).Data.EEG.raw.value);    %filtering EEG data
end
figure
plot(freq,abs(fft(feeg(1,1:length(freq)))).^2);
title('PSD of eeg')
%% Extract delta waves
w1=0.5/(0.5*fs);  %delta waves lower cut off frequency
w2=3/(0.5*fs); %delta waves lower cut off frequency
wn=[w1,w2]; %setting up filter coefficients
[balpha,aalpha]=butter(2,wn); %creating filter
for i=1:size(feeg,1)
    deeg(i,:)=filtfilt(balpha,aalpha,feeg(i,:));    %getting delta waves
end
%% Getting the frequency response of delta waves
freq=0:(fs/length(eeg(1).Data.EEG.raw.value)):40; 
for i=1:size(feeg,1)
    deegdft(i,:)=fft(deeg(i,:));  %applying fft on delta waves
end
%% Extract theta waves
w1=3/(0.5*fs);  %theta waves lower cut off frequency
w2=8/(0.5*fs); %theta waves lower cut off frequency
wn=[w1,w2]; %setting up filter coefficients
[balpha,aalpha]=butter(2,wn); %creating filter
for i=1:size(feeg,1)
    teeg(i,:)=filtfilt(balpha,aalpha,feeg(i,:));    %getting theta waves
end
%% Getting the frequency response of theta waves
freq=0:(fs/length(eeg(1).Data.EEG.raw.value)):40; 
for i=1:size(feeg,1)
    teegdft(i,:)=fft(teeg(i,:));  %applying fft on theta waves
end
%% Extract alpha waves
w1=8/(0.5*fs);  %alpha waves lower cut off frequency
w2=13/(0.5*fs); %alpha waves lower cut off frequency
wn=[w1,w2]; %setting up filter coefficients
[balpha,aalpha]=butter(2,wn); %creating filter
for i=1:size(feeg,1)
    aeeg(i,:)=filtfilt(balpha,aalpha,feeg(i,:));    %getting alpha waves
end
%% Visualise alpha waves in time domain
% figure
% plot(aeeg(1,:));
% title('Alpha EEG');
% xlabel('Time');
% ylabel('Amplitude');
%% Getting the frequency response of delta waves
freq=0:(fs/length(eeg(1).Data.EEG.raw.value)):40; 
for i=1:size(feeg,1)
    aeegdft(i,:)=fft(aeeg(i,:));  %applying fft on alpha waves
end
%% Extract beta waves
w1=13/(0.5*fs);  %theta waves lower cut off frequency
w2=38/(0.5*fs); %theta waves lower cut off frequency
wn=[w1,w2]; %setting up filter coefficients
[balpha,aalpha]=butter(2,wn); %creating filter
for i=1:size(feeg,1)
    beeg(i,:)=filtfilt(balpha,aalpha,feeg(i,:));    %getting theta waves
end
%% Getting the frequency response of beta waves
freq=0:(fs/length(eeg(1).Data.EEG.raw.value)):40; 
for i=1:size(feeg,1)
    beegdft(i,:)=fft(beeg(i,:));  %applying fft on beta waves
end
%% Extract gamma waves
w1=38/(0.5*fs);  %theta waves lower cut off frequency
w2=42/(0.5*fs); %theta waves lower cut off frequency
wn=[w1,w2]; %setting up filter coefficients
[balpha,aalpha]=butter(2,wn); %creating filter
for i=1:size(feeg,1)
    geeg(i,:)=filtfilt(balpha,aalpha,feeg(i,:));    %getting theta waves
end
%% Getting the frequency response of beta waves 
for i=1:size(feeg,1)
    geegdft(i,:)=fft(geeg(i,:));  %applying fft on beta waves
end
%% Clear variables not required anymore
clear a b balpha aalpha i q wn w1 w2
%% Visualise power spectral density of eeg components
for i=1:size(eeg,2)
    response(:,:,i)=[deegdft(i,:); teegdft(i,:); aeegdft(i,:); beegdft(i,:); geegdft(i,:)];
    signal(:,:,i)=[deeg(i,:); teeg(i,:); aeeg(i,:); beeg(i,:); geeg(i,:)];
end
figure
hold on
plot(freq,abs(response(1,1:length(freq),2)).^2); title('PSD of components of  EEG');
plot(freq,abs(response(2,1:length(freq),2)).^2);
plot(freq,abs(response(3,1:length(freq),2)).^2);
plot(freq,abs(response(4,1:length(freq),2)).^2);
plot(freq,abs(response(5,1:length(freq),2)).^2);
legend('Delta','Theta','Alpha','Beta','Gamma')
% change response's index accordingly to choose the wave you want to see
xlabel('Frequency');
ylabel('Power');
%% Extracting statistical features from extracted 