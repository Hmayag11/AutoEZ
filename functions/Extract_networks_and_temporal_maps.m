function [networks_en,networks_bk,temporal_map]=Extract_networks_and_temporal_maps(P,ids0,k)
% [networks_en,networks_bk,temporal_map]=Extract_networks_and_temporal_maps(P_sb,ids0,k)
% Extracts the networks and the temporal maps of 30 repetitions using NNMF
%
% INPUTS: 
%       P: DMD spectral features
%       ids0: ids of the modes to process
%       k: number of componetns to extract (we used k=2 in the manuscript)
%
% OUTPUTS:
%       networks_en: epileptogenic network
%       networks_bk: background network
%       temporal_map: temporal maps

%% Repeating the nnmf analysis 30 times and collect them in an array

temporal_map=[];
networks_en=[];
networks_bk=[];
for i=1:30
    [m,r,w]=size(P);
    % average powers of all frequency components per band
    Phi_current=abs(squeeze(mean(P(:,ids0,:),2)));

    % Extract the repeating spatial configurations and their temporal
    % coefficients using NNMF
    [nets,coeffs,normbest]=nnmf(abs(Phi_current),k);
    nets=nets';

    [val,ind1]=max(coeffs);% max is used to find the initial temoral map
    c1=sum(ind1==1);% frequency of spatial configuration 1
    c2=sum(ind1==2);% frequency of spatial configuration 2
    if c1<c2 % reorder so that the less frequent is the second row
        nets_new(1,:)=nets(2,:);
        nets_new(2,:)=nets(1,:);
        ind1=abs(ind1-3);
    else
        nets_new(1,:)=nets(1,:);
        nets_new(2,:)=nets(2,:);
    end
    index=kmeans(coeffs',2);% kmeans is applied to smooth the temporal maps
    c3=norm(ind1'-index);
    c4=norm(ind1'-(3-index));
    if c4<c3 % correct the indices of the kmeans to match those ordered by the nnmf max indices
        index=abs(index-3);
    end
    temporal_map=[temporal_map;ind1];
    networks_bk=[networks_bk;nets_new(1,:)];
    networks_en=[networks_en;nets_new(2,:)];
end
    
% c1,c2,c3,c4 help 
% 1. detemrine the less frequent network (epileptogenic)
% 2. correct the indices of the kmeans 
% 3.always set the first network to bakcground and the second one to epileptogenic