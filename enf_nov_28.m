% %enf.m function file
% % 
% fs = 44100;  
% blocksize = fs*16; 
% zero_pad =0; 
% overlap = 0.5; 
% w = hanning( blocksize,'periodic' ); 
% user_freq = 240; 
% %Test File: recording.wav 
% [recordingdata, ~] = audioread('ground truth 2.wav');  
% y = recordingdata; 

%%****** Start of the function ******** %%
function [y1, y2] = enf_nov_28(y,fs,blocksize,zero_pad,overlap,w,user_freq)
indent = 1; %indent can also be referred to as row number
product = round(overlap*blocksize);  
y = y.';
x = y; %input data  
%%%%%%%%%%%%% Hamming Window %%%%%%%%%%%%%%%%%% 
%w = hann(blocksize,'periodic'); 
w = w.';  
% %%%%%%%%%%initial block%%%%%%%%%%%%%%%%%%%%%%%
 xm(indent,:) = [x(1:blocksize)];   
 hm(indent,:) = xm(indent,:).*w;
 zpm(indent,:) = horzcat(hm(1,:),zeros(1,zero_pad));
 start_indent = blocksize - (product - 1); 
 end_indent = start_indent + (blocksize-1); 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%Subsequant blocks after %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% loop goes off of the initial block and adds on it to
additional_blocks = floor(length(x)/product) - 1; %estimating the
%amount of rows/loops that will be needed   
%additional_blocks = 27;
%% 
for i = 1:additional_blocks
    indent = indent + 1;  %increment indent/row number 
if (end_indent > length(x))  % If statement for adding zeroes if the input data cannot evenly divide into windows
    diff = end_indent - length(x); %figuring how many zeroes to add at the end if empty space
    end_indent = length(x);   
    xcond = [x(start_indent:end_indent)]; 
    xm(indent,:) = [xcond zeros(1,diff)];  
    break;
else  
   xm(indent,:) = [x(start_indent:end_indent)]; %main line that creates the blocks 
   
end  
%Hanning windowing of segmented input.
hm(indent,:) = xm(indent,:).*w;

%%%%%%%%%%%%%%%%%5%%%%%%%%%Zero Padding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i==1)
    zpm(indent-1,:) = horzcat(hm(i,:),(zeros(1,zero_pad)));
else
    zpm(i, :) = horzcat(hm(i,:),(zeros(1,zero_pad)));
end

start_indent = end_indent - (product - 1); %Look at the end value of the previous row, skip to the product indent,  
% take that as the beginning for the next row
end_indent = start_indent + (blocksize - 1); %start from the start value and go up to size blocksize, minusing one  
%since working with indents 
end 

%%%%%%%%%%%%%%%%%% Zero Padding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Zero padded value of 
zpm(additional_blocks, :) = horzcat(hm(i,:),(zeros(1,zero_pad)));  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Short Time Fourier Time %%%%%%%%%%%%%%%%%%%%% 
fftm = abs(fft(zpm,[],2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%% User Defined Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
user_f1 = user_freq - 1; 
user_f2 = user_freq + 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
newsize = size(zpm()); 
N = newsize(1,2) ;   %updated block size after adding zeroes
k_indent = ( round(((user_f1*N)/fs)+1) : 1 : round(((user_f2*N)/fs)+1) );  
ENF_region = fftm(:,k_indent);      
surf(ENF_region); 
%%example: E = (k,m) or (m,k) where k => 33 , m => 27
%%m => block index, aka row number 

%%k => (16*2) + 1  
%%k => ((F(k) * N)/Fs) + 1
%%F(k) => F(k1) = F(kenf)minus1 && F(k2) = F(kenf)plus1 (bassically
%%119 and 121, user_f1 and user_f2)
%%fs => sampling frequency 
%%F(#) = ((k-1)Fs)/N , # is the indent 
%%F(k) already created (user_f1 and user_f2) 
%%Therefore E(k,m) => ENF region (k,m)


%%weighted_energy    

%(N/Fs * 2) + 1, if N = 705600 and Fs = 44100 
% ((705600 / 44100) * 2 ) + 1 => 33 <- the k value in E(k,m) 
% m -> bins/frames -> 27 => specific to recording, got it do be 28
% by variable indent 
k = ((N/fs)*2)  + 1; %supposed to be 33
F_k = linspace(user_f1, user_f2, k); %after calc should be 28 rows
                                     %and 33 colums
% Taking the F_k now need to repeat it for all of the rows that
% ENF_region has to ensure that both matrices are the same size and
% the F_k is being summed with each row of the ENF_region.
% Therefore will need to use repmat  
row_size = size(ENF_region); 
indent = row_size(1,1);
F_k_rep = repmat(F_k, indent , 1); %%where indent => the row size of the
                     %ENF_region (28)
                     % need to limit to 1 copy per row or else
                     % won't be 33   
                     % currently making copies a 1 x 33 array for
                     % indent rows (indent = 28)
a = F_k_rep .* ENF_region; 
%format of sum => sum(matrix,dim)
weighted_energy = sum(a,2) ./ sum(ENF_region,2);



%%max energy
[~,max_energy] = max(ENF_region,[],2); 
% y1 => maximum energy 
y1 = ( max_energy .* (fs/N) ) + user_f1;
% y2 => weighted energy
y2 = weighted_energy;
