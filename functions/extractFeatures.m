function [P,freq_mean,L]=extractFeatures(ieeg,perc,Fs,r_input,offset,annotation)
% [P,freq_mean,LL]=extractFeatures(ieeg,perc,Fs,r_input,offset,spike_annotation)
% Extracts features using the Dynamic Mode Decomposition using a sliding
% window approach
%
% INPUTS: 
%       ieeg: channels x time samples
%       perc: percentage overlap
%       Fs: sampling frequency
%       r_input: number of DMD modes
%       offset: controls the time in seconds before and after a point in time. it
%       helps specify the size fo the window
%       annotation: the time points with spike or ripple annotations
%
% OUTPUTS:
%       P: DMD spectral features
%       freq_mean: mean frequency of modes
%       L: time-window labels, contianing spikes/ripples or not
% 

%initialize variables as empty
window=[];
Phi_all=[];
F=[];
L=[];
FF=[];
PP=[];
P=[]';

% extract DMD spectra using a sliding window approach

[n,m]=size(ieeg);
loc=400:perc*Fs:m-450;% get rid of the boundaries
st1=offset*Fs;% offset seconds before
end1=offset*Fs;% offset seconds after
for jj=1:length(loc)
    disp(['Processing window',' ',num2str(jj)])
    
    window=[(ieeg(:,round(loc(jj))- round(st1):round(loc(jj))+round(end1)))];% choose a window
    h=ceil(2*size(window,2)/size(window,1))+1;% number of time delay embeddings
    
    for r=r_input
        [Phi, omega, lambda, b, Xdmd, time_dynamics,V_r,U_r,W_r,Atilde,S_r,X1,X2] = DMDfull(window, [],h,r,1); % apply DMD on current time-window
        freq=abs(omega)/(2*pi)*Fs;
    end

    Phi_all(:,:,jj)=Phi(1:n,:); % contains the spatial modes
    F(:,jj)=freq; % contians the frequency of oscillatory eigenfunctions

    % label each window whether it contains spikes or not, based on
    % annotatons
    V = round(annotation*Fs);
    N = round(loc(jj))- round(st1):round(loc(jj))+round(end1);
    for i=1:length(V)
        [c2(i) index(i)] = min(abs(N-V(i)));
    end
    if min(c2)==0
        L(jj)=1;
    else
        L(jj)=0;
    end
end

% postprocess by ordering Phi's and frequencies in increasing order
[~,~,p]=size(Phi_all);
for i=1:p
    PP=[abs(squeeze(Phi_all(:,:,i))) ;F(:,i)'];
    B = sortrows(PP',n+1);
    FF(:,i)=B(:,n+1);
    P(:,:,i)=B(:,1:n)';
end

freq_mean=mean(FF')'; %mean frequency across windows
end