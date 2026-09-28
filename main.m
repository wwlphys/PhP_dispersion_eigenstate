fidin=fopen('model.input','r');
itype=fscanf(fidin,'%d',1);
if itype==1
    disp('calculate bulk PhP dispersion')
    calculate_dispersion
elseif itype==2
    disp('calculate slab PhP dispersion')
    calculate_surface_dispersion
else
    disp('Error: not supported itype')
    itype
end