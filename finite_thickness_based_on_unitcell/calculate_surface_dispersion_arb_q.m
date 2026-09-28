clear
fidin=fopen('iq.txt','r');
iq=fscanf(fidin,'%d',1)
fclose(fidin);
firsttime=1;
if firsttime==0
    dm_preprocess;
    phonopy_yaml_preprocess;
    OUTCAR_preprocess;
else
    disp('Skip preprocess!');
end
read_dm;
%test_eig;
read_VnM;
read_Z;
read_FORCE;
%expand_FORCE;
%test_dm;
%findrank2
%find0
%comp3factorwp
LargeDet_minw_arb_q
exit
