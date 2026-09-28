clear
disp('processing qpoints.yaml');
dm_preprocess;
disp('reading clean_qpoints.txt');
read_dm;
%test_eig;
disp('processing phonopy.yaml');
phonopy_yaml_preprocess;
disp('reading clean_phonopy.txt');
read_VnM;
OUTCAR_preprocess;
read_Z;
findzeropointrealallq_insert  % 计算所有q点，用于一维色散关系或二维等频线
%findzeropointRANDq_insert  % 随机顺序地计算qpoints.yaml中的q点，用于多个节点共同计算三维等频面
