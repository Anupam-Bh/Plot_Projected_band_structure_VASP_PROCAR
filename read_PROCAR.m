%% read PROCAR (single PROCAR needed: if splitted make a combined PROCAR)
clc
clear all
tic
filename='PROCAR';
fid=fopen(filename,'r');
% read #kpts, #bands and #ions
while feof(fid)==0
    S=fgetl(fid);
    str=sprintf(S);
    clear X S
    X=strsplit(str,{' ','\t'});
    if strcmp(X(1),'#')==1
        nkpt=str2double(X(4))
        nband=str2double(X(8))
        nion=str2double(X(12))
        break;
    end
end
clear str;
fclose(fid);

fid=fopen(filename,'r');
ikpt=0;
iband=0;
iion=0;
while feof(fid)==0
    S=fgetl(fid);
    str=sprintf(S);
    clear X S
    X=strsplit(str,{' '});
    if length(X)>=2 && strcmp(X(2),'k-point')==1
        ikpt=ikpt+1;
        kpt(ikpt,:)=str2double(X(5:7));
        iband=0;
        if mod(ikpt,10)==0
         fprintf(' %d kpoints done: ',ikpt);
         toc
        end
    end
    if strcmp(X(1),'band')==1
        iband=iband+1;
        iion=0;
        energy(ikpt,iband)=str2double(X(5));
        occ(ikpt,iband)=str2double(X(8));
    end
    if length(X)>=2 && iion<nion && str2double(X(2))-1==iion && strcmp(X(1),'band')==0
        iion=iion+1;
        spd(ikpt,iband,iion,:)=str2double(X(3:length(X)-1));
    end
end
fclose(fid);
dim=size(spd);
norb=dim(4);

param=[nkpt,nband,nion,norb];
save('MATLAB_param','param');
save('MATLAB_spd','spd');
save('MATLAB_eigs','energy');
save('MATLAB_occupations','occ');
save('MATLAB_kpoints','kpt');
