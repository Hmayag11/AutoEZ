%% add the required paths
addpath(genpath('sample_data'))  % add path of the data
addpath(genpath('functions'))   %add path of funcitons
addpath(genpath('anatomy'))  %add path of funcitons

% clear and close all
clear
close all

% set figflag to 1 to pop-up result figures, otherwise set it to 0
figflag=1;

%% load the data
load('sample_data_P31.mat')

%% plot the MRI in 3D and overlay resection and electrodes
V = niftiread('MRI_P31.nii');
V = rot90(V);
[X,Y,Z] = meshgrid(-0.1445:0.0011705:0.1558-0.0011705,-0.1363:0.00103:0.1276-0.00103,-0.1415:0.00115:0.1534-0.00115);
if figflag
    figure
    xslice =0;
    yslice = 0;
    zslice = 0;
    colormap(bone)
    % plot the central slices in the three directions
    h=slice(X,Y,Z,double(V),xslice,yslice,zslice);
    view(180,0)
    set(h,'edgecolor','none')
    % plots the resection volume
    hold on
    bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
    trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',1,'Edgecolor',[0 1 0]);
    % plots the channels
    hold on
    [X1,Y1,Z1] = sphere(100);
    val=ones(size(channel_coordinates,1),1)/700;
    for i=1:length(val)
        surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','w','Edgecolor','w','FaceAlpha',0.5)
    end
    axis square
    title('MRI, Resection Volume, and Electrodes of Patient 31')
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    set(gca,'zticklabel',{[]})
    set(gcf,'color','w');
    h = zeros(2, 1);
    h(1) = scatter(NaN,NaN,'og','filled');
    h(2) = scatter(NaN,NaN,'ok');
    l=legend(h, '\color{green} Resection','Electrodes');
    fontsize(l,14,'points')
    set(gcf, 'Position', get(0, 'Screensize'));
    xlim([-0.1445,0.1539])
    zlim([-0.14,0.15])
    saveas(gcf, './results/MRIandResectionandElectrodes_P31.png')
end

%% plot the data segment
 channels=cell(size(channel_coordinates,1));
tickMarks = {'YTick',[]};
[n,m]=size(ieeg);
if figflag
    figure;
    for i=1:n
        plot(time,ieeg(i,:)-0.0015*i,'k','linewidth',1.5);
        hold on
    end
    hold on
    xline(spike_annotation,'r')
    xlim([0 50])
    ylim([-0.0015*i 0])
    set(gcf,'color','w')
    yticks(fliplr(-0.0015*(1:n)))
    for i=1:n
        channels{i}=num2str(n-i+1);
    end
    xlabel('time (s)')
    ylabel('channel ID')
    yticklabels(channels)
    set(gca,'fontsize',8)
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Sample Data for Patient 31')
    set(gcf, 'Position', get(0, 'Screensize'));
end

%% spike band analysis
%filter the data in sb band [1-80]Hz
ieeg_filtered_spike=zeros(size(ieeg));
f1=1;
f2=80;
[b,a]=butter(4,[2*f1/(Fs),2*f2/(Fs)]);
for i=1:n
    ieeg_filtered_spike(i,:)=filtfilt(b,a,double(ieeg(i,:)));
end
% plot the filtered signal
if figflag
    figure;
    for i=1:n
        plot(time,ieeg_filtered_spike(i,:)-0.0015*i,'k','linewidth',1.5);
        hold on
    end
    hold on
    xline(spike_annotation,'r')
    xlim([0 50])
    ylim([-0.0015*i 0])
    set(gcf,'color','w')
    yticks(fliplr(-0.0015*(1:n)))
    for i=1:n
        channels{i}=num2str(n-i+1);
    end
    xlabel('time (s)')
    ylabel('channel ID')
    yticklabels(channels)
    set(gca,'fontsize',8)
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Sample segment in the spike band for Patient 31')
    set(gcf, 'Position', get(0, 'Screensize'));
end
%% Using the sliding window approach,extract the features
perc=0.05; % percentage overlap
r_input=50;% number of modes
offset=0.125; % determines how much (in seconds) before and after a point in time (determines the window size, here 0.25 seconds total)
[P_sb,freq_mean,Lsb]=extractFeatures(ieeg_filtered_spike,perc,Fs,r_input,offset,spike_annotation);

%% identify bands of modes, collect the IDs using the mean frequency of all windows
ids0=getIndices(freq_mean,r_input);

%% extract the networks in the six subbands of sb
k=2; % we assumed two networks and set k=2
networks_new=cell(7);
temporal_map_new=zeros(7,length(Lsb));
for bd=1:6
    ids00=ids0{bd};
    [networks_en,networks_bk,temporal_map]=Extract_networks_and_temporal_maps(P_sb,ids00,k);
    % averaging to find consistent networks and temporal maps
    networks_new{bd}=[mean(networks_bk);mean(networks_en)];
    temporal_map_new(bd,:)=mean(temporal_map)>1.5;
end

%% ripple band analysis (same flow as spike)
% filter the data in sb band [1-80]Hz
ieeg_filtered_ripple=zeros(size(ieeg));
f1=80;
f2=250;
[b,a]=butter(4,[2*f1/(Fs),2*f2/(Fs)]);
for i=1:n
    ieeg_filtered_ripple(i,:)=filtfilt(b,a,double(ieeg(i,:)));
end
% plot the data in the ripple band
if figflag
    figure;
    for i=1:n
        plot(time,ieeg_filtered_ripple(i,:)-0.0004151*i,'k','linewidth',1.5);
        hold on
    end
    xlabel('time (s)')
    hold on
    xline(ripple_annotation,'r')
    xlim([0 50])
    set(gcf,'color','w')
    yticks(fliplr(-0.0004151*(1:n)))
    for i=1:n
        channels{i}=num2str(n-i+1);
    end
    xlabel('time (s)')
    ylabel('channel ID')
    yticklabels(channels)
    a = get(gca,'XTickLabel');
    set(gca,'fontsize',8)
    title('Sample Data in ripple band')
    set(gcf, 'Position', get(0, 'Screensize'));
end
%% Using the sliding window approach,extract the features
perc=0.05; %percentage overlap
r_input=100;% number of modes
offset=0.075;% determines how much (in seconds) before and after a point in time (determines the window size, here 0.15 seconds total)
[P_rb,freq_mean,Lrb]=extractFeatures(ieeg_filtered_ripple,perc,Fs,r_input,offset,ripple_annotation);

%% Identify networks and temporal maps using NNMF
ids00=1:100;% all modes are used for ripple band
[networks_en,networks_bk,temporal_map]=Extract_networks_and_temporal_maps(P_rb,ids00,k);
% averaging to find consistent networks and temporal maps
networks_new{7}=[mean(networks_bk);mean(networks_en)];
temporal_map_new(7,:)=mean(temporal_map)>1.5;

%% plot temporal maps
end1=50*Fs;
if figflag
    figure;
    subplot(4,1,1)
    for i=11:20  % For illustration purposes, we used 10 channels with spikes/ripples
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
    imagesc(temporal_map_new(1:6,:));colormap(CustomColormap)%1:19.5*end1/Fs
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
    imagesc(temporal_map_new(7,:));colormap(CustomColormap)
    set(gca,tickMarks{:});
    yticks(1)
    yticklabels({'\it rb'})
    title('Temporal map in \itrb')
    xlabel('window ID')
    sgtitle('Patient 31 temporal maps of 50s sample iEEG data')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, './results/TemporalMap50sec_P31.png')
end

%% plot temporal maps for the first 10 seconds
st1=1;
end1=10*Fs;% % 10 controls the number of seconds to plot
if figflag
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
    yticks(1)
    yticklabels({'\it rb'})
    title('Temporal map in \itrb')
    xlabel('window ID')
    sgtitle('Patient 31 temporal maps of 10s sample iEEG data')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, './results/TemporalMap10sec_P31.png')
end

%% plot the epileptogenic and background networks
% choose band to analyze:  1= delta; 2=theta;3=alpha;4=beta;5=gamma;6=sb;7=rb%
bd_current=2; 

if figflag
    figure
    subplot(1,2,1)
    % plots the three slices of the MRI
    xslice =0;
    yslice = 0;
    zslice = 0;
    colormap(bone)
    h=slice(X,Y,Z,double(V),xslice,yslice,zslice);
    view(180,0)
    set(h,'edgecolor','none')
    % plots the resection volume
    hold on
    bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
    trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',1,'Edgecolor',[0 1 0]);
    % plots the network with spheres whose radii are proportional to the power
    [X1,Y1,Z1] = sphere(100);
    hold on
    val=ones(size(channel_coordinates,1),1)/300;
    network=networks_new{bd_current};
    val=val.*network(2,:)';% 1 for background 2 for epileptogenic
    for i=1:length(val)
        surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','r','Edgecolor','r','FaceAlpha',0.5)
    end
    axis square
    xlim([-0.1445,0.1539])
    zlim([-0.14,0.15])
    title('\color{red} Epileptogenic Network','FontSize', 24)
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    set(gca,'zticklabel',{[]})
  
    % same process to plot the bakcground network in another subplot
    subplot(1,2,2)
    % plots the three slices of the MRI
    xslice =0;
    yslice = 0;
    zslice = 0;
    colormap(bone)
    h=slice(X,Y,Z,double(V),xslice,yslice,zslice);
    view(180,0)
    set(h,'edgecolor','none')
    % plots the resection volume
    hold on
    bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
    trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',1,'Edgecolor',[0 1 0]);
    %plots the network with spheres whose radii are proportional to the power
    [X1,Y1,Z1] = sphere(100);
    val=ones(size(channel_coordinates,1),1)/300;
    network=networks_new{bd_current};% 2 for theta
    val=val.*network(1,:)';% 1 for background 2 for epileptogenic
    for i=1:length(val)
        surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','b','Edgecolor','b','FaceAlpha',0.5)
    end
    axis square
    xlim([-0.1445,0.1539])
    zlim([-0.14,0.15])
    title('\color{blue} Background Network','FontSize', 24)
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    set(gca,'zticklabel',{[]})
    set(gcf,'color','w');
    sgtitle('Patient 31 Networks')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, './results/Networks_P31.png')
end

%% Overlayed networks
if figflag
    figure
    xslice =0;
    yslice = 0;
    zslice = 0;
    colormap(bone)
    % plots the three slices of the MRI
    h=slice(X,Y,Z,double(V),xslice,yslice,zslice);
    view(180,0)
    set(h,'edgecolor','none')
    % plots the resection volume
    hold on
    bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
    trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',1,'Edgecolor',[0 1 0]);
    %plots the network with spheres whose radii are proportional to the power
    [X1,Y1,Z1] = sphere(100);
    hold on
    val=ones(size(channel_coordinates,1),1)/300;
    network=networks_new{bd_current};
    val=val.*network(2,:)';% 1 for background 2 for epileptogenic
    for i=1:length(val)
        surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','r','Edgecolor','r','FaceAlpha',0.5)
    end
    axis square
    xlim([-0.1445,0.1539])
    zlim([-0.14,0.15])
    title('Epileptogenic Network')
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    set(gca,'zticklabel',{[]})

    % perform the same process to plot the background network superimposed
    hold on
    xslice =0;
    yslice = 0;
    zslice = 0;
    colormap(bone)
    h=slice(X,Y,Z,double(V),xslice,yslice,zslice);
    view(180,0)
    set(h,'edgecolor','none')
    hold on
    bound=boundary(resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),0.7);
    trisurf(bound,resection_coordinates(:,1),resection_coordinates(:,2),resection_coordinates(:,3),'Facecolor',[ 0 1 0],'FaceAlpha',1,'Edgecolor',[0 1 0])
    [X1,Y1,Z1] = sphere(100);
    val=ones(size(channel_coordinates,1),1)/300; % scaling to fit coordinates dimensions
    network=networks_new{bd_current};% 2 for theta
    val=val.*network(1,:)';%% 1 for background 2 for epileptogenic
    for i=1:length(val)
        surf(X1*val(i)+channel_coordinates(i,1),Y1*val(i)+channel_coordinates(i,2),Z1*val(i)+channel_coordinates(i,3),'Facecolor','b','Edgecolor','b','FaceAlpha',0.5)
    end
    axis square
    xlim([-0.1445,0.1539])
    zlim([-0.14,0.15])
    title('Epileptogenic and Backlground Networks of Patient 31')
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    set(gca,'zticklabel',{[]})
    set(gcf,'color','w');
    hold on;
    h = zeros(2, 1);
    h(1) = scatter(NaN,NaN,'or','filled');
    h(2) = scatter(NaN,NaN,'ob','filled');
    l=legend(h, '\color{red} Epileptogenic','\color{blue} Background');
    fontsize(l,14,'points')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, './results/Networks Overlayed_P31.png')
end


%% properties of the networks
% computes the focality, overlap, and distance from resection for both
% networks
Fnet=zeros(2,7);
Ores=zeros(2,7);
Dres=zeros(2,7);
res_thr=10;
for bd=1:7
    network=networks_new{bd};
    for j=1:2 % j =1:2 for background and epileptogenic
        U=network(j,:)>mean(network(j,:))+std(network(j,:));
        Fnet(j,bd)=1000*mean(pdist(channel_coordinates(U>0,:)));
        coord=channel_coordinates(U>0,:);
        dist=distiance_electrodes_to_resection(coord,resection_coordinates);
        dist_fin=double(dist(:,4));
        res_elec_map=dist_fin<res_thr;
        Ores(j,bd)=100*sum(res_elec_map)/length(res_elec_map); % use it for outcome confusion matrix
        Dres(j,bd)=mean(dist_fin);
    end
end
disp(Fnet)
disp(Ores)
disp(Dres)
%% AUCSOZ
% computes AUC of each network to predict the SOZ
AUCSOZ=zeros(2,7);
for bd=1:7
    network=networks_new{bd};
    for j=1:2
        U=network(j,:).*(network(j,:)>0.3);
        [~,~,~,AUCSOZ(j,bd)] = perfcurve(soz_channels,U,1);
    end
end
disp(AUCSOZ)

%% AUCRES
% computes AUC of each network to predict the resection
AUCRES=zeros(2,7);
for bd=1:7
    network=networks_new{bd};
    for j=1:2 
        coord=channel_coordinates;
        dist=distiance_electrodes_to_resection(coord,resection_coordinates);
        dist_fin=double(dist(:,4));
        res_elec_map=dist_fin<res_thr;
        U=network(j,:).*(network(j,:)>0.3);
        [~,~,~,AUCRES(j,bd)] = perfcurve(res_elec_map,U,1);
    end
end
disp(AUCRES)

%% AUC-IED
% computes AUC of each network's activity in the temporal map to detect the
% IEDs/ripple(for rb)
AUCIED=zeros(2,7);
temporal_map=temporal_map_new;
for bd=1:7
    map=temporal_map(bd,:);
    map_bk=(map==0);
    map_epil=(map==1);
    if bd<7
        [~,~,~,AUCIED(1,bd)] = perfcurve(Lsb,0.9*map_bk,1);
        [~,~,~,AUCIED(2,bd)] = perfcurve(Lsb,0.9*map_epil,1);
    else
        [~,~,~,AUCIED(1,bd)] = perfcurve(Lrb,0.9*map_bk,1);
        [~,~,~,AUCIED(2,bd)] = perfcurve(Lrb,0.9*map_epil,1);
    end
end
disp(AUCIED)

%% summary results
summary=[Fnet;Ores;Dres;AUCSOZ;AUCRES;AUCIED];

%% Plot Network Properties
x=[1 2];
legends={'EN','BN'};
color= ['r','g'];
if figflag
    figure
    subplot(2,3,1)
    b=bar(x,Fnet(:,bd_current)');
    b.FaceColor = 'flat';
    b.CData(1,:) =[0    0   1];
    b.CData(2,:) = [1    0   0];
    title('Fnet (mm)')
    ylim([0 80])
    text(1:length(x),Fnet(:,bd_current)',num2str(round(Fnet(:,bd_current))),'vert','bottom','horiz','center');

    subplot(2,3,2)
    b=bar(Ores(:,bd_current));
    b.FaceColor = 'flat';
    b.CData(1,:) =[0    0   1];
    b.CData(2,:) = [1    0   0];
    title('Ores (%)')
    ylim([0 120])
    text(1:length(x),Ores(:,bd_current)',num2str(round(Ores(:,bd_current))),'vert','bottom','horiz','center');

    subplot(2,3,3)
    b=bar( Dres(:,bd_current));
    title('Dres (mm)')
    b.FaceColor = 'flat';
    b.CData(1,:) =[0    0   1];
    b.CData(2,:) = [1    0   0];
    ylim([0 50])
    text(1:length(x),Dres(:,bd_current)',num2str(round(Dres(:,bd_current))),'vert','bottom','horiz','center');

    subplot(2,3,4)
    b=bar( AUCRES(:,bd_current));
    b.FaceColor = 'flat';
    b.CData(1,:) =[0    0   1];
    b.CData(2,:) = [1    0   0];
    title('AUC-RES')
    ylim([0 1.2])
    text(1:length(x), AUCRES(:,bd_current)',num2str(round( AUCRES(:,bd_current),2)),'vert','bottom','horiz','center');

    subplot(2,3,5)
    b=bar( AUCSOZ(:,bd_current));
    b.FaceColor = 'flat';
    b.CData(1,:) =[0    0   1];
    b.CData(2,:) = [1    0   0];
    title('AUC-SOZ')
    ylim([0 1.2])
    text(1:length(x), AUCSOZ(:,bd_current)',num2str(round( AUCSOZ(:,bd_current),2)),'vert','bottom','horiz','center');

    subplot(2,3,6)
    b=bar(x,AUCIED(:,bd_current)');
    b.FaceColor = 'flat';
    b.CData(1,:) =[0    0   1];
    b.CData(2,:) = [1    0   0];
    ylim([0 1.2])
    title('AUC-IED')
    text(1:length(x), AUCIED(:,bd_current)',num2str(round(AUCIED(:,bd_current),2)),'vert','bottom','horiz','center');
    set(gcf,'color','w')
    sgtitle('Patient 31 Network Properties in the \theta band ({\color{blue}BN\color{black} & \color{red}EN}) ')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, './results/NetworkProperties_P31.png')
end


