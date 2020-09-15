%ENFproject.m   

%format of function
% [y1, y2] = enf(y,fs,blocksize,zero_pad,overlap,w,user_freq) 

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 60; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 60Hz of Ground 1'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 60Hz of recording 1'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 60Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 60Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 60Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 60Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 60Hz of Ground 1'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 60Hz of recording 1'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 60Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 60Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%120 HZ with Audio Ground 1 and Recording 1

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 120; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 120Hz of Ground 1'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 120Hz of recording 1'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 120Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 120Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 120Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 120Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 120Hz of Ground 1'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 120Hz of recording 1'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 120Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 120Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%180Hz with audio 1

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 180; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 60Hz of Ground 1'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 60Hz of recording 1'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 180Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 180Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 180Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 180Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 180Hz of Ground 1'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 180Hz of recording 1'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 180Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 180Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%240 Hz with Audio recording 1 and ground 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%180Hz with audio 1

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 240; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 240Hz of Ground 1'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 240Hz of recording 1'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 240Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 240Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 240Hz of Ground 1'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 240Hz of recording 1'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 240Hz of Ground 1'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 240Hz of recording 1'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 240Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 240Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 60; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth 2.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording 2.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 60Hz of Ground 2'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 60Hz of recording 2'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 60Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 60Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 60Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 60Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 60Hz of Ground 2'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 60Hz of recording 2'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 60Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 60Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%120 HZ with Audio Ground 2 and Recording 2

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 120; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth 2.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording 2.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 120Hz of Ground 2'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 120Hz of recording 2'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 120Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 120Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 120Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 120Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 120Hz of Ground 2'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 120Hz of recording 2'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 120Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 120Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%180Hz with audio 2

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 180; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth 2.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording 2.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 60Hz of Ground 2'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 60Hz of recording 2'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 180Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 180Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 180Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 180Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 180Hz of Ground 2'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 180Hz of recording 2'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 180Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 180Hz'); 
xlabel('Lags');
ylabel('Correlation');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%240 Hz with Audio recording 1 and ground 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%240Hz with audio 2

clear;

Fs = 44100;  
BlockSize = Fs*16; 
Zeropad = 0; 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequency = 240; 

%Reference File: ground truth.wav 
[gndtruthdata, ~] = audioread('ground truth 2.wav'); 

%Test File: recording.wav 
[recordingdata, ~] = audioread('recording 2.wav');  

%Max and Weighted energy of ground truth.wav   
figure(); 
subplot(2,1,1);
[y1gndmax, y2gndweighted] = enf_nov_28(gndtruthdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency);
title('ENF Surface Plot of Ground Wave No Preprocess At 240Hz of Ground 1'); 
%Max and Weighted energy of recording.wav  
subplot(2,1,2);
[y1recmax, y2recweighted] = enf_nov_28(recordingdata,Fs,BlockSize,Zeropad,Overlap,Window,Frequency); 
title('ENF Surface Plot of Recording Wave No Preprocess At 240Hz of recording 1'); 


%Transposing the data  
%y2recweighted = y2recweighted.';
%y2gndweighted = y2gndweighted.';

%% Normalized Cross corelation
%Normalized cross correlation (max energy) 
% a[m] and b[m] are the arrays 
% l is the lag factor 
%  ua and ub are the mean values of a[m] and b[m] 

%to ensure that both arrays can be used in sum, both must be the same size 
% xcorr solves this issue by finding the lagging vector and adding zeros to
% it thus making it same size 

% if x and y is a multidimensional array, xcorr returns autocorr
% and cross-corr as columns in a matrix 
% lag indices are returned as a vector, indices where lag occurs?
% r is the cross-correlation/autocorrelation vector  

%Max mean
ua_gndmax = mean(y1gndmax); % will return a row vector of mean values for each column 
ub_recmax = mean(y1recmax);  


gndmax_test = transpose(y1gndmax);
recmax_test = transpose(y1recmax);


[rmax,lagmax] = xcorr(y1gndmax, y1recmax);    


%Plot Max energy
figure();
subplot(2,1,1);
plot (gndmax_test);
title(' Ground max at 240Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recmax_test);
title(' Recorded max at 240Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  

%Weighted mean  
ua_gndweight = mean(y2gndweighted); % will return a row vector of mean values for each column 
ub_recweight = mean(y2recweighted);

gndweight_test = transpose(y2gndweighted);
recweight_test = transpose(y2recweighted);

%Plot weighted energy
figure();
subplot(2,1,1);
plot (gndweight_test);
title(' Ground weighted enery at 240Hz of Ground 2'); 
xlabel('Lags');
ylabel('Amplitudes');  
subplot(2,1,2);
plot (recweight_test);
title(' Recorded weighted enery at 240Hz of recording 2'); 
xlabel('Lags');
ylabel('Amplitudes');  



% Figure is plotted in the decimated section for more easier comparison
%figure();
[rweight,lagweight] = xcorr(y2gndweighted - ua_gndweight, y2recweighted - ub_recweight);   
subplot(2,1,2);
%plot (lagweight, rweight);
%title('Correalation vs. Lag, gnd and rec weight w/ NO preproc'); 
%xlabel('Lags');
%ylabel('Correlation');   



%% Normalized corelation with pre-processing 

%need to decimate sound files by 100  
% need to filter then downsample further 

% Before filter must down sample the signal 
% passband of 150 Hz
% stopband of 500 Hz 
% FDA Tool -> IIR filter -> stop 2 pass in ppt? -> export to matlab
% -> use filtfilt command 
% -> make sure to pick the one that uses SOS and G tool   

load('lowpassfilter.mat'); 


gndtruth_filtered = filtfilt(SOS,G,gndtruthdata);
recordingdata_filtered = filtfilt(SOS,G,recordingdata);

% now need to decimate/downsample the filtered function by 100
% use decimate function y = decimate (x,r) 
% decimates input signal x by factor of r 

gndtruthfilt_decimate = decimate (gndtruth_filtered, 100); 
recordingdatafilt_decimate = decimate (recordingdata_filtered, 100);

Fsp = Fs/100;  
BlockSize = Fsp*16; 
Zeropad = 16384 - (16*Fsp); 
Overlap = 0.5; 
Window = hanning( BlockSize,'periodic' ); 
Frequencyp = Frequency;  


%Max and Weighted energy of preprocessed ground truth.wav  
figure(); 
subplot(2,1,1);
[y1gndmax_preproc, y2gndweighted_preproc] = enf_nov_28(gndtruthfilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp);
title('ENF Surface Plot of Ground Wave Preprocessed At 240Hz of Ground 2'); 

%Max and Weighted energy of preprocessed recording.wav 
subplot(2,1,2);
[y1recmax_preproc, y2recweighted_preproc] = enf_nov_28(recordingdatafilt_decimate,Fsp,BlockSize,Zeropad,Overlap,Window,Frequencyp); 
title('ENF Surface Plot of Recording Wave Preprocessed At 240Hz of recording 2'); 



%Weighted preprocessed mean  
ua_gndweight_preproc = mean(y2gndweighted_preproc); % will return a row vector of mean values for each column 
ub_recweight_preproc = mean(y2recweighted_preproc);  



figure(); 
[rweight_preproc,lagweight_preproc] = xcorr(y2gndweighted_preproc - ua_gndweight_preproc, y2recweighted_preproc - ub_recweight_preproc);   
subplot(2,1,1); 
plot (lagweight, rweight);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ NO preproc At 240Hz'); 
xlabel('Lags');
ylabel('Correlation');   

subplot(2,1,2); 
plot (lagweight_preproc, rweight_preproc);
title('Correalation vs. Lag, gnd2 and rec2 weight w/ preproc At 240Hz'); 
xlabel('Lags');
ylabel('Correlation');   

