clear
path(path,'/public2/home/openclaw/public_code/finite_thickness_based_on_unitcell/eigenvector');
path(path,'/public2/home/openclaw/public_code/finite_thickness_based_on_unitcell');
path(path,'/public2/home/openclaw/public_code');
firsttime=0;
if firsttime==0
    dm_preprocess;
    phonopy_yaml_preprocess;
    OUTCAR_preprocess;
    band_yaml_preprocess;
else
    disp('Skip preprocess!');
end
read_dm;
%test_eig;
read_VnM;
read_Z;
read_FORCE;
read_displace;
%expand_FORCE;
%test_dm;
%findrank2
%find0
%comp3factorwp
LargeDet_minw_eigenvector_ps
%exit
