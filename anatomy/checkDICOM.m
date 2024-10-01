%% this script loads DICOM files, info, and plots the DICOM in 3D. 
% Load DICOM
filename='MRI_P1.dcm';
X = dicomread(filename);

% Plot a slice,
s=150;%slice number
figure;imagesc(X(:,:,1,150))

% Extract DICOM info
info = dicominfo( filename);
info

% Plot DICOM in 3D
sx = 1;
sy = 1;
sz = 2.5;
A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
tform = affinetform3d(A);
J = permute(X, [1 2 4 3]);
vol = volshow(squeeze(J(:,:,:,1)),renderingStyle="SlicePlanes",Transformation=tform);