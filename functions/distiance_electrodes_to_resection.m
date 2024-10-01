function [file_output] =distiance_electrodes_to_resection(channel_coordinates,resection_coordinates)

% [dipole_file_output] =compare_dipole_to_iEEG_Location(dipole_file,electrode_coord,distance_coverage)
% Compute the distance between different sets of coordinates (we mainly use
% this to find the distances of the electrodes from the resection volume)
%
% INPUTS: 
%       channel_coordinates: set of electrodes (channels)
%       resection_coordinates: set of electrodes (resection volume)
%       both Nx3 matrix (double), where N is the number of electrodes (columns are x-y-z coordinates)
%
% OUTPUTS:
%       file_output: the coordinates fo the channels and their distance from resection in mm


if isa(channel_coordinates,'struct')
    channel_coordinates=struct2dataset(channel_coordinates');
    for i=1:size(channel_coordinates,1)
        loc=(channel_coordinates.Loc{i});
        loc_x(i,1)=loc(1);
        loc_y(i,1)=loc(2);
        loc_z(i,1)=loc(3);
    end
    channel_coordinates.loc_x=loc_x;
    channel_coordinates.loc_y=loc_y;
    channel_coordinates.loc_z=loc_z;
end
if isa(channel_coordinates,'double')
    channel_coordinates=dataset(channel_coordinates(:,1),channel_coordinates(:,2),channel_coordinates(:,3));
    channel_coordinates.Properties.VarNames{1}='loc_x';
    channel_coordinates.Properties.VarNames{2}='loc_y';
    channel_coordinates.Properties.VarNames{3}='loc_z';
end
if isa(resection_coordinates,'dataset')
    electrode_coord_new(:,1)=resection_coordinates.loc_x;
    electrode_coord_new(:,2)=resection_coordinates.loc_y;
    electrode_coord_new(:,3)=resection_coordinates.loc_z;
    resection_coordinates=electrode_coord_new;
end


loc_x = channel_coordinates.loc_x;
loc_y = channel_coordinates.loc_y;
loc_z = channel_coordinates.loc_z;

loc=[loc_x loc_y loc_z];
n_all=size(loc,1);


min_dist_iEEG=[]; min_dist_elec=[];

n_all=size(loc,1);
for i=1:n_all
    dip_dist=[];
    for j=1:size(resection_coordinates,1)
        dip_dist(j,1)=pdist([loc(i,:); resection_coordinates(j,:)]);
    end
    [min_d, min_ind]=min(dip_dist);
    min_dist_iEEG(i,1)=min_d;
    min_dist_elec(i,1)=min_ind;
end
distance_iEEG=min_dist_iEEG;

file_output=channel_coordinates;
file_output.distance_from_resection=distance_iEEG*1000;



end