iq=1;
qvec=[qp(iq,1) qp(iq,2) qp(iq,3)];
%qvec=[0 0*vb(2,2) 0.5*vb(3,3)]
%DMq
%w3=eig(DM1);
%real(sqrt(w3))/(2*pi);
fid=fopen('testdm0.txt','w');
for i=0:0
    qvec=[vb(1,1)*0 vb(2,2)*0 vb(3,3)*0];
    DMq
    %disp(DM1(1,4:6));
    %disp(DM1(48,27));
    clear w3
    w3=eig(DM1);
    freq=real(sqrt(w3))/(2*pi)
    for j=1:nband
        plot(i,freq(j),'*');
        hold on
    end
    for j=1:nband
        for k=1:nband
            fprintf(fid,'%g  %g  ',real(DM1(j,k)),imag(DM1(j,k)));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);