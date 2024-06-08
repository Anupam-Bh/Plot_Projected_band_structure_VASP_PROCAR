%% Plot spd decomposed band-structure
% comment different sections as per needed
clear all
close all
clc
load('MATLAB_eigs.mat');
load('MATLAB_kpoints.mat');
load('MATLAB_param.mat');
load('MATLAB_occupations.mat');
load('MATLAB_spd.mat');
nkpt=param(1)
nband=param(2)
nion=param(3)
norb=param(4)
A=fileread('OUTCAR');
B=strfind(A,'E-fermi');
efermi=str2double(A(B+10:B+19))
% efermi=input('enter fermi level: ');
Xpoints=1:nkpt;
%% Enter high symmetry points

vertices_xval=[1 30  60 90 120 150 180 210  240 270 300 330 360]
vertices_char={'\Gamma' 'X' 'S' 'Y' '\Gamma' 'Z' 'U' 'R' 'T' 'Z|X' 'U|Y' 'T|S' 'R'}
% vertices_xval=[1 1000 2000]
% vertices_char={'L''' '\Gamma' 'L'};
%  vertices_xval=[1 100 200 300 400 500 600 700 800 900 1000]
%  vertices_char={'\Gamma','L','B1|B','Z','\Gamma','X|Q','F','P1','Z','L','P'};  

%% plot bandstructure
figure(1)
for i=1:nband
    plot(1:nkpt,energy(:,i)-efermi,'k','LineWidth',1);
    hold on
end
set(gca,'xlim',[0 max(Xpoints)],'ylim',[-7 7],'Xtick',vertices_xval,'Xticklabel',vertices_char,'Xgrid','on','Ygrid','on',...
       'Fontweight','normal','Fontsize',18,'Fontname','arial');
ylabel('E-E_{f} (eV)');
% title('Band-structure');
box on;
pbaspect([1.5 1 1]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
%% plot atomic contribution
atom_proj_3d=sum(spd,4);
% Types of atoms : Read POSCAR
filename='POSCAR';
fid=fopen(filename,'r');
count=0;
while feof(fid)==0
    S=fgetl(fid);
    count=count+1;
    if count == 6
        str=sprintf(S);
        atom_name=strsplit(str,{' ','\t'});
        atom_name=atom_name(2:end)
    elseif count == 7
        str=sprintf(S);
        typeatoms=strsplit(str,{' ','\t'});
        typeatoms=typeatoms(2:end)
        break;
    end
end
atom_proj_type=zeros(nkpt,nband,length(typeatoms));
for nn=1:length(typeatoms)
    atom_proj_type(:,:,nn)=sum(atom_proj_3d(:,:,(nn-1)*cell2sym(typeatoms(nn))+1:nn*cell2sym(typeatoms(nn))),3);
end

% with maximum atomic contributions
atom_proj=zeros(nkpt,nband);
atom_proj1=zeros(nkpt,nband);
for i=1:nkpt
    for j=1:nband
        [xx,yyy]=max(atom_proj_type(i,j,:));
        atom_proj(i,j)=yyy;atom_proj1(i,j)=xx;
    end
end
figure(2)
for i=1:nband
      scatter(1:nkpt,energy(:,i)-efermi,atom_proj1(:,i).*5,atom_proj(:,i),'filled')
%     scatter(1:nkpt,energy(:,i)-efermi,5,atom_proj(:,i),'filled')
hold on
end
set(gca,'xlim',[0 max(Xpoints)],'ylim',[-15 7],'Xtick',vertices_xval,'Xticklabel',vertices_char,'Xgrid','on','Ygrid','on',...
       'Fontweight','normal','Fontsize',18,'Fontname','times');
ylabel('E-E_{f} (eV)');
colormap(jet);
title('Atomic contribution(max contribution)');
c=colorbar;
c.Ticks = [1:length(typeatoms)];
c.TickLabels = atom_name;
box on;
pbaspect([1.5 1 1]);

        
%% Plot spd projected BS (Can also be used for atomic contribution) 
% give the # of atoms in form of an array e.g. [2 4 7]
atoms=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16] 
% give the # of orbitals in form of an array e.g. [3 5 16]
% s=[1]
% py, pz and px =[2 3 4]
% dxy  dyz  dz2  dxz  dx2-y2 = [5 6 7 8 9]
% fy3x2  fxyz  fyz2  fz3  fxz2  fzx2  fx3=[10 11 12 13 14 15 16]
lmproj=[1 2 3 4 5 6 7 8 9]
% lmproj=input('Enter which orbitals are to be projected: ');
%summing atoms which are to be projected
at=repmat({':'},1,4);
at{3}=atoms;
% at{4}=lmproj;
proj=spd(at{:});
proj_spd=squeeze(sum(proj,3));
totproj=1:16; totproj(lmproj)=[];
proj_spd(:,:,totproj)=0;
% determining spdf character
proj_2d=zeros(nkpt,nband);
for i=1:nkpt
    for j=1:nband
        %[maxi,in]=max(proj_spd(i,j,lmproj));
        sproj(i,j)=proj_spd(i,j,1);
        pproj(i,j)=sum(proj_spd(i,j,2:4));
        dproj(i,j)=sum(proj_spd(i,j,5:9));
        fproj(i,j)=sum(proj_spd(i,j,10:16));
        allproj(i,j)=sum(proj_spd(i,j,1:16));
        
        [maxi,in]=max([sproj(i,j) pproj(i,j) dproj(i,j) fproj(i,j)]);
         proj_2d_weight(i,j)=maxi;
        %index=lmproj(in);
        if maxi==0
            proj_2d(i,j)=0;
        elseif in==1
            proj_2d(i,j)=5;   %s orbital
        elseif in==2
            proj_2d(i,j)=10;  %p orbital
        elseif in==3
            proj_2d(i,j)=15;  %d orbital
        elseif in==4
            proj_2d(i,j)=20;  %f orbital
        end
    end
end

figure(3)
for i=1:nband
% scatter(1:nkpt,energy(:,i)-efermi,10,proj_2d(:,i),'filled')
scatter(1:nkpt,energy(:,i)-efermi,5*(proj_2d_weight(:,i))+1e-3,proj_2d(:,i),'filled')
% scatter(1:nkpt,energy(:,i)-efermi,10*(allproj(:,i)),[0.6 .5 .7],'filled')
% scatter(1:nkpt,energy(:,i)-efermi,10*(sproj(:,i))+1e-4,'red','filled')
% hold on
% scatter(1:nkpt,energy(:,i)-efermi,10*(pproj(:,i))+1e-4,'green','filled')
% hold on
% scatter(1:nkpt,energy(:,i)-efermi,10*(dproj(:,i))+1e-4,'blue','filled')
% hold on
% scatter(1:nkpt,energy(:,i)-efermi,50*(fproj(:,i)),'orange','filled')
hold on
end
% vertices_char={'L''' '\Gamma' 'L'}
set(gca,'xlim',[0 max(Xpoints)],'ylim',[-15 7],'Xtick',vertices_xval,'Xticklabel',vertices_char,'Xgrid','on','Ygrid','on',...
       'Fontweight','normal','Fontsize',18,'Fontname','arial');
ylabel('E-E_{f} (eV)');
colormap(jet);
% text(1500,0.1,'s-type','Color',[0 0.5 0.4],'FontSize',22)
% text(1500,-0.05,'p-type','Color',[0.64 0.8 0.38],'FontSize',22)
% text(1500,-0.14,'d-type','Color',[1.0 0.12 0],'FontSize',14)
dd=colorbar;
dd.Ticks = [5 10 15 20];
dd.TickLabels = {'s' 'p' 'd' 'f'};
% caxis([5 12])
title('Orbital projected band-structure');
box on;
pbaspect([1.5 1 1])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
%% plot occupations
% figure(4)
% for i=1:nband
%     scatter(1:nkpt,energy(:,i)-efermi,5,occ(:,i),'filled')
%     hold on
% end
% set(gca,'xlim',[0 max(Xpoints)],'ylim',[-1 1],'Xtick',vertices_xval,'Xticklabel',vertices_char,'Xgrid','on','Ygrid','on',...
%        'Fontweight','bold','Fontsize',17,'Fontname','times');
% ylabel('E-E_{f} (eV)');
% colormap(jet);
% box on;
% title('occupation');
% caxis('auto');

