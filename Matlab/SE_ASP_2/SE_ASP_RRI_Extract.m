clc;clear;
%{
%Vijay data
load('C:\Users\Vijay\Google Drive\Year 4\Spring\ASP\Matlab\ECG\RAW_Vijay.mat')

%RRI data extraction%
a0=300;
a1=2.45e5;

b0=2.9e5;
b1=4.9e5;

c0=5.2e5;
c1=7.43e5;


trial_1 = data(a0:a1);
trial_2 = data(b0:b1);
trial_3 = data(c0:c1);


%[RRI_trial_1,fs_RRI] = ECG_to_RRI(trial_1,fs);   %ampthresh=0.1
%[RRI_trial_2,fs_RRI] = ECG_to_RRI(trial_2,fs);   %ampthresh=0.7
%[RRI_trial_3,fs_RRI] = ECG_to_RRI(trial_3,fs);   %ampthresh=0.75
%}

%{
%Will Smith data
load('C:\Users\Vijay\Google Drive\Year 4\Spring\ASP\Matlab\ECG\RAW_WillSmith2.mat')

%RRI data extraction%
a0=250;
a1=2.42e5;

b0=2.8e5;
b1=4.85e5;

c0=5.23e5;
c1=7.35e5;

trial_1 = data(a0:a1);
trial_2 = data(b0:b1);
trial_3 = data(c0:c1);

%[RRI_trial_1,fs_RRI] = ECG_to_RRI(trial_1,fs);   %ampthresh=0.2
%[RRI_trial_2,fs_RRI] = ECG_to_RRI(trial_2,fs);   %ampthresh=0.2
%[RRI_trial_3,fs_RRI] = ECG_to_RRI(trial_3,fs);   %ampthresh=0.14
%}
