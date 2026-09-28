clear
dm_preprocess;
read_dm;    % qpoints.yaml
%test_eig;
phonopy_yaml_preprocess;
read_VnM;        % phonopy.yaml
OUTCAR_preprocess;
read_Z;          % OUTCAR.txt
if 1  %  要不要跟声子本征态比较
    band_yaml_preprocess;
    read_displace  % band.yaml
else
    ph_dspl=zeros(nband,nband);
    disp('Not comparing with phonon. ovlp.txt is set to zero.');
end
findeigenvector3D_largeDI;
