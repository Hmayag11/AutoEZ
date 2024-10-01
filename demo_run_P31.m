%% add the required paths
addpath('C:\Users\hpartamian-temp\CBD_Team Dropbox\PHDUTA PHDUTA\Behnam\Toolboxes\brainstorm\toolbox\core\')
addpath('C:\Users\hpartamian-temp\Desktop\forJournalFinal\sample_data')
addpath('C:\Users\hpartamian-temp\Desktop\forJournalFinal\functions')

clear all
close all
tic
%% load the data
load('data_pt2.mat')

%% open the MRI in 3D, resection volume and channels
% you may need to add the following brainstorm  addpath('C:\Users\~~~~~~~\Toolboxes\brainstorm\toolbox\core\')
openfig('MRI_Pt2.fig') 
view(1.7,-0.5)

%% click anywhere in the 3D figure then run this. 
% If the MRI is loaded from brainstorm, no need to click

% plots the resections volume
hold on
bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',0.44,'Edgecolor',[0 1 0])

% plots the channels
hold on
[X1,Y1,Z1] = sphere(100);
[m,n]=size(channel_coordinates);
val=ones(size(channel_coordinates,1),1)/1000;
for i=1:length(val)
surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','k','Edgecolor','k','FaceAlpha',0.5)
end

%% plot segment of data
tickMarks = {'YTick',[]};

[n,m]=size(ieeg);
figure;
for i=1:n
    plot(time,ieeg(i,:)-0.00151*i,'k','linewidth',1.5);
    hold on
end
hold on
xline(spike_annotation,'r')
xlim([0 50])
ylim([-0.00151*n 0])
set(gcf,'color','w')
% set(gca,tickMarks{:});
yticks(fliplr(-0.00151*(1:n)))
for i=1:n
channels{i}=num2str(n-i+1);
end

yticklabels(channels)
a = get(gca,'XTickLabel');
set(gca,'fontsize',8)
title('Sample Data')
%% spike band analysis
%filter the data in sb band [1-80]Hz
f1=1;
f2=80;

[b,a]=butter(4,[2*f1/(Fs),2*f2/(Fs)]);
for i=1:n
    ieeg_filtered_spike(i,:)=filtfilt(b,a,double(ieeg(i,:)));
end

figure;
for i=1:n
    plot(time,ieeg_filtered_spike(i,:)-0.0019151*i,'k','linewidth',1.5);
    hold on
end
%
xlabel('time (s)')
hold on
xline(spike_annotation,'r')
xlim([0 50])
% ylim([-0.125 0])
set(gcf,'color','w')
yticks(fliplr(-0.0019151*(1:n)))

for i=1:n
channels{i}=num2str(n-i+1);
end

yticklabels(channels)
a = get(gca,'XTickLabel');
set(gca,'fontsize',8)
title('Sample Data in spike band')

%% Using the sliding window approach,extract the features 
perc=0.1; % percentage overlap 
r_input=50;% number of modes
offset=0.125; % determines how much (in seconds) before and after a point in time (determines the window size, here 0.25 seconds total)
[P_sb,freq_mean,Lsb]=extractFeatures(ieeg_filtered_spike,perc,Fs,r_input,offset,spike_annotation);

%% identify bands of modes, collect the IDs using the mean frequency of all windows 
ids0=getIndices(freq_mean,r_input);

%%
k=2;
networks_new=[];
networks_new=[];
temporal_map_new=[];

for bd=1:6
 ids00=ids0{bd};
[networks_en,networks_bk,temporal_map]=Extract_networks_and_temporal_maps(P_sb,ids00,k);
% averaging to find consistent networks and temporal maps
networks_new{bd}=[mean(networks_bk);mean(networks_en)];
temporal_map_new=[temporal_map_new;mean(temporal_map)];
end

%% ripple band analysis (same flow as spike)
%filter the data in sb band [1-80]Hz
f1=80;
f2=250;
[b,a]=butter(4,[2*f1/(Fs),2*f2/(Fs)]);  
 for i=1:n
    ieeg_filtered_ripple(i,:)=filtfilt(b,a,double(ieeg(i,:)));
end

 figure;
for i=1:n
    plot(time,ieeg_filtered_ripple(i,:)-0.0004151*i,'k','linewidth',1.5);
    hold on
end
%
xlabel('time (s)')
hold on
xline(ripple_annotation,'r')
xlim([0 50])
% ylim([-0.125 0])
set(gcf,'color','w')
 yticks(fliplr(-0.0004151*(1:n)))
for i=1:n
channels{i}=num2str(n-i+1);
end
yticklabels(channels)
a = get(gca,'XTickLabel');
set(gca,'fontsize',8)
title('Sample Data in ripple band')

%% Using the sliding window approach,extract the features 
perc=0.1; %percentage overlap 
r_input=100;% number of modes
offset=0.075;% determines how much (in seconds) before and after a point in time (determines the window size, here 0.15 seconds total)
[P_rb,freq_mean,Lrb]=extractFeatures(ieeg_filtered_ripple,perc,Fs,r_input,offset,ripple_annotation);

%% Identify networks and temporal maps using NNMF
ids00=1:100;% All modes are used for ripple band
[networks_en,networks_bk,temporal_map]=Extract_networks_and_temporal_maps(P_rb,ids00,k);
networks_new{7}=[mean(networks_bk);mean(networks_en)];
temporal_map_new=[temporal_map_new;mean(temporal_map)];

%% plot temporal maps
end1=50*Fs;
figure;
subplot(4,1,1)
for i=11:20
    plot(time(1:end1),ieeg_filtered_spike(i,1:end1)-0.00151*i,'k','linewidth',1.5);
    hold on
end
hold on
xline(spike_annotation,'r')
xlim([0 50])
set(gcf,'color','w')
set(gca,tickMarks{:});
xlabel('time (s)')
ylabel('Channels')
title('Sample data in spike band with spike annotations(\itred)')

subplot(4,1,2)
imagesc(temporal_map_new(1:6,:)>1.5);colormap(CustomColormap)%1:19.5*end1/Fs
set(gca,tickMarks{:});
yticks([1 2 3 4 5 6])
yticklabels({'\delta','\theta','\alpha' , '\beta' ,'\gamma', '\it sb'})
title('Temporal map in \it{sb} subbands')
xlabel('window ID')

subplot(4,1,3)
for i=11:20
    plot(time(1:end1),ieeg_filtered_ripple(i,1:end1)-0.0004151*i,'k','linewidth',1.5);
    hold on
end
hold on
xline(spike_annotation,'r')
xlim([0 50])
set(gcf,'color','w')
set(gca,tickMarks{:});
xlabel('time (s)')
ylabel('Channels')
title('Sample data in spike band with ripple annotations(\itred)')

subplot(4,1,4)
imagesc(temporal_map_new(7,:)>1.5);colormap(CustomColormap)
set(gca,tickMarks{:});
yticks([1])
yticklabels({'\it rb'})
title('Temporal map in \itrb')
xlabel('window ID')
%% plot temporal maps for the first 10 seconds
 
st1=1;
end1=10*Fs;
figure;
subplot(4,1,1)
for i=11:20
    plot(time(st1:end1),ieeg_filtered_spike(i,st1:end1)-0.00151*i,'k','linewidth',1.5);
    hold on
end
hold on
xline(spike_annotation,'r')
xlim([0.13 10])
set(gcf,'color','w')
set(gca,tickMarks{:});
xlabel('time (s)')
ylabel('Channels')
title('Sample data in spike band with spike annotations(\itred)')

subplot(4,1,2)
imagesc(temporal_map_new(1:6,1:196));colormap(CustomColormap)%1:19.5*end1/Fs
set(gca,tickMarks{:});
yticks([1 2 3 4 5 6])
yticklabels({'\delta','\theta','\alpha' , '\beta' ,'\gamma', '\it sb'})
title('Temporal map in \it{sb} subbands')
xlabel('window ID')

subplot(4,1,3)
for i=11:20
    plot(time(st1:end1),ieeg_filtered_ripple(i,st1:end1)-0.0004151*i,'k','linewidth',1.5);
    hold on
end
hold on
xline(spike_annotation,'r')
xlim([0.13 10])
set(gcf,'color','w')
set(gca,tickMarks{:});
xlabel('time (s)')
ylabel('Channels')
title('Sample data in spike band with ripple annotations(\itred)')

subplot(4,1,4)
imagesc(temporal_map_new(7,1:196));colormap(CustomColormap)
set(gca,tickMarks{:});
yticks([1])
yticklabels({'\it rb'})
title('Temporal map in \itrb')
xlabel('window ID')
% subplot(3,1,1)
% for i=1:7
% plot(networks{i}')
% hold on;
% end
% plot(soz_channels,'k','linewidth',2)

%% plotting the networks

f=openfig('MRI_Pt2.fig')
view(1.7,-0.5)
 
%% click anywhere in the 3D figure then run this to plot the epileptogenic network
bd=2; %2= theta 3=alpha

hold on
bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',0.44,'Edgecolor',[0 1 0])
%
[X1,Y1,Z1] = sphere(100);
hold on


val=ones(size(channel_coordinates,1),1)/300;
network=networks_new{bd};% 2 for theta
val=val.*network(2,:)';%% 1 for background 2 for epileptogenic
for i=1:length(val)
surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','r','Edgecolor','r','FaceAlpha',0.5)
end
%%

f=openfig('MRI_Pt2.fig')
view(1.7,-0.5)
% axes
%% click anywhere in the 3D figure then run this to plot the background network
bd=2; %2= theta 3=alpha

hold on
bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',0.44,'Edgecolor',[0 1 0])

[X1,Y1,Z1] = sphere(100);
[m,n]=size(channel_coordinates);
val=ones(size(channel_coordinates,1),1)/300;
network=networks_new{bd};% 2 for theta
val=val.*network(1,:)';%% 1 for background 2 for epileptogenic
for i=1:length(val)
surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','b','Edgecolor','b','FaceAlpha',0.5)
end


%% properties of the networks
addpath('C:\Users\hpartamian-temp\OneDrive - UT Arlington\Microsoft Teams Chat Files\')
res_thr=10;% 10mm around resection
for bd=1:7;
    network=networks_new{bd};%
    for j=1:2 %background and epil
        U=network(j,:)>0.2;%mean(network(j,:))+std(network(j,:))
        Fnet(j,bd)=1000*mean(pdist(channel_coordinates(U>0,:)));
        coord=channel_coordinates(U>0,:);
        dist=distiance_electrodes_to_resection(coord,resection_coordinates);
        dist_fin=double(dist(:,4));
        res_elec_map=dist_fin<res_thr;
        Ores(j,bd)=100*sum(res_elec_map)/length(res_elec_map); % use it for outcome confusion matrix
        Dres(j,bd)=mean(dist_fin);
    end
end
Fnet
Ores
Dres
%% AUCSOZ

for bd=1:7;
    network=networks_new{bd};%
    for j=1:2 %background and epil
    U=network(j,:).*(network(j,:)>0.2);%%mean(network(j,:))+std(network(j,:)))
    [X,Y,t1,AUCSOZ(j,bd)] = perfcurve(soz_channels,U,1);
    end
end
AUCSOZ

%% AUCRES
for bd=1:7;
    network=networks_new{bd};%
    
    for j=1:2 %background and epil
        U=network(j,:)>0.3;%%resultsPost{1, 1}.thresholds(pt,bd);
        coord=channel_coordinates;
        dist=distiance_electrodes_to_resection(coord,resection_coordinates);
        dist_fin=double(dist(:,4));
        res_elec_map=dist_fin<res_thr;
    U=network(j,:).*(network(j,:)>0.3);
    [X,Y,t1,AUCRE(j,bd)] = perfcurve(res_elec_map,U,1);
    end
end
AUCRE

%% AUC-IED
temporal_map=temporal_map_new;
for bd=1:7
    map=temporal_map(bd,:);%
    map_bk=(map==0);
    map_epil=(map==1);
    if bd<7
    [X,Y,t1,AUC_IED(1,bd)] = perfcurve(Lsb,0.9*map_bk,1);
     [X,Y,t1,AUC_IED(2,bd)] = perfcurve(Lsb,0.9*map_epil,1);
    else
     [X,Y,t1,AUC_IED(1,bd)] = perfcurve(Lrb,0.9*map_bk,1);
     [X,Y,t1,AUC_IED(2,bd)] = perfcurve(Lrb,0.9*map_epil,1);
    end
end
AUC_IED

%%
summary=[Fnet;Ores;Dres;AUCSOZ;AUCRE;AUC_IED]

%% Plot Network Properties 

bd=2 %% choose band  2=theta
x=[1 2];
legends={'EN','BN'};

color= ['r','g'];
figure
subplot(2,3,1)
b=bar(x,Fnet(:,bd)');
b.FaceColor = 'flat';
b.CData(1,:) =[0    0   1];
b.CData(2,:) = [1    0   0];
title('Fnet (mm)')
ylim([0 80])
text(1:length(x),Fnet(:,bd)',num2str(round(Fnet(:,bd))),'vert','bottom','horiz','center'); 

subplot(2,3,2)
b=bar( Ores(:,bd))
b.FaceColor = 'flat';
b.CData(1,:) =[0    0   1];
b.CData(2,:) = [1    0   0];
title('Ores (%)')
ylim([0 120])
text(1:length(x),Ores(:,bd)',num2str(round(Ores(:,bd))),'vert','bottom','horiz','center'); 

subplot(2,3,3)
b=bar( Dres(:,bd));
title('Dres (mm)')
b.FaceColor = 'flat';
b.CData(1,:) =[0    0   1];
b.CData(2,:) = [1    0   0];
ylim([0 50])
text(1:length(x),Dres(:,bd)',num2str(round(Dres(:,bd))),'vert','bottom','horiz','center'); 

subplot(2,3,4)
b=bar( AUCRE(:,bd));
b.FaceColor = 'flat';
b.CData(1,:) =[0    0   1];
b.CData(2,:) = [1    0   0];
title('AUC-RES')
ylim([0 1.2])
text(1:length(x), AUCRE(:,bd)',num2str(round( AUCRE(:,bd),2)),'vert','bottom','horiz','center'); 

subplot(2,3,5)
b=bar( AUCSOZ(:,bd))
b.FaceColor = 'flat';
b.CData(1,:) =[0    0   1];
b.CData(2,:) = [1    0   0];
title('AUC-SOZ')
ylim([0 1.2])
text(1:length(x), AUCSOZ(:,bd)',num2str(round( AUCSOZ(:,bd),2)),'vert','bottom','horiz','center'); 

subplot(2,3,6)
b=bar(x,AUC_IED(:,bd)');
b.FaceColor = 'flat';
b.CData(1,:) =[0    0   1];
b.CData(2,:) = [1    0   0];
ylim([0 1.2])
title('AUC-IED')
text(1:length(x), AUC_IED(:,bd)',num2str(round(AUC_IED(:,bd),2)),'vert','bottom','horiz','center'); 
set(gcf,'color','w')
sgtitle('Network Properties in the \theta band ({\color{blue}BN\color{black} & \color{red}EN})')
% sgtitle(['\fontsize{16}black {\color{magenta}magenta \color[rgb]{0 .5 .5}teal \color{red}red} black again'])

toc