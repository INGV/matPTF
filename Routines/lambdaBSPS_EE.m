function [lambda_BS,lambda_PS,lambda_PS_zone,tetra_vol]=lambdaBSPS_EE(hypo_utm,utmzone_hypo,cov,SD,Npts_vs_Ntetra,Nvref_vref_ratio,moho_pts,Nz,mesh_subd,out_tetra_vol)
disp('----- running: lambdaBSPS_EE ------')

%% parameters %%



%% ---- EarlyEst Hypo and Cov setting in utm -----------------%
hypox=hypo_utm(:,1);
hypoy=hypo_utm(:,2);
hypoz=hypo_utm(:,3);

S=eig(cov);
l1=SD*sqrt(S(1,1)); % major semi-axis in km
l2=SD*sqrt(S(2,1)); % intermediate semi-axis in km
l3=SD*sqrt(S(3,1)); % minor semi-axis in km
%--------------------------------------------%
VolEll=(4/3)*pi*l1*l2*l3; % ellipsoid theoretical volume

if (license('test', 'curve_fitting_toolbox'))
    Npts_Mw=ceil(Npts_vs_Ntetra(VolEll.*Nvref_vref_ratio));
else
    Npts_Mw=ceil(Npts_vs_Ntetra(1).*(VolEll.*Nvref_vref_ratio).^Npts_vs_Ntetra(2));
end

Npts_Mw=max(10,Npts_Mw);
disp(['--> lambdaBSPS_EE:  Npts_Mw = ' num2str(Npts_Mw)]);


%% ---- Ellipsoid + Tetrahedra transform. -----%
[xp,yp,zp] = plot_gaussian_ellipsoid_noplot([hypox hypoy -hypoz*1000],cov,SD,Npts_Mw); close all;
x2=xp(:);
y2=yp(:);
z2=zp(:);
[vertices,ia,ic]=unique([x2,y2,z2],'rows');
vertices_ll=vertices;
[vertices_ll(:,2),vertices_ll(:,1)]=utm2ll(vertices(:,1),vertices(:,2),utmzone_hypo);
tetra = delaunayn(vertices_ll); % Generate delaunay triangulization
%--------------------------------------------%


%% ---- Check the Ellipsoid Volume ------------%
Ntetra=size(tetra,1);
volume=zeros(Ntetra,1);
for i=1:size(tetra,1)
    volume(i,1)=abs(det([[vertices(tetra(i,:),1)';vertices(tetra(i,:),2)';vertices(tetra(i,:),3)'];1 1 1 1]))/6;
end


Vol_tetra=sum(volume); % ellipsoid volume from tetrahedra

Vol_diff_perc=(VolEll-Vol_tetra)./VolEll*100; % Percentage difference between theoretical and tetrahedra volumes
%--------------------------------------------%


%% ---- Tetrahedra Baricenters ----------------%
baric_tetra=zeros(Ntetra,3);

baric_tetra(:,1)=mean(reshape(vertices_ll(tetra,1),Ntetra,4),2);
baric_tetra(:,2)=mean(reshape(vertices_ll(tetra,2),Ntetra,4),2);
baric_tetra(:,3)=mean(reshape(vertices_ll(tetra,3),Ntetra,4),2);

pts=baric_tetra;
pts(:,3)=-pts(:,3)/1000;
[ptsx,ptsy]=ll2utm(pts(:,2),pts(:,1),utmzone_hypo);


if (out_tetra_vol==1)
    tetra_vol=cell(2,1);
    tetra_vol{1,1}=baric_tetra;
    tetra_vol{2,1}=volume;
else
    tetra_vol=[];
end

%--------------------------------------------%

%--------------Moho grid init------------------%

depth_moho_pts=griddata(moho_pts(1).moho_par(:,1),moho_pts(1).moho_par(:,2),moho_pts(1).moho_par(:,3),pts(:,1),pts(:,2));

%---- Subduction Zones Meshes ---------------%

dist_pts_mesh_eu=zeros(Ntetra,1);

iPS_sep=cell(Nz,1);
iPS_eff_sep=cell(Nz,1);
distances_sep=cell(Nz,1);

for i=1:Nz
        
    %---- Euclidean Distance between the Ellipsoid points and the Subd Mesh triangles baricenters
    Ntri=size(mesh_subd(4).mesh_par{i,1},1);
    tri_mesh=zeros(Ntri,3);
    
    tri_mesh(:,1)=mean(reshape(mesh_subd(5).mesh_par{i,1}(mesh_subd(4).mesh_par{i,1}(1:Ntri,2:end),1),Ntri,3),2);
    tri_mesh(:,2)=mean(reshape(mesh_subd(6).mesh_par{i,1}(mesh_subd(4).mesh_par{i,1}(1:Ntri,2:end),1),Ntri,3),2);
    tri_mesh(:,3)=mean(reshape(mesh_subd(3).mesh_par{i,1}(mesh_subd(4).mesh_par{i,1}(1:Ntri,2:end),4),Ntri,3),2);
    
    
    [dist_pts_mesh_eu, iid]=pdist2([tri_mesh(:,1),tri_mesh(:,2),-tri_mesh(:,3)],...
        [ptsx(:,1),ptsy(:,1),pts(:,3)*1000],'euclidean','smallest',1);
    dist_pts_mesh_eu=dist_pts_mesh_eu'/1000;
    distances_sep{i,1}=dist_pts_mesh_eu;
    
    %---- Find PS points ------------------------%
    tsh=10; % buffer (in km) around the subduction interface
    iPS_sep{i,1}=find(distances_sep{i,1}<=tsh);
    iPS_eff_sep{i,1}=find(-pts(:,3)<0 & distances_sep{i,1}<=tsh);
    %--------------------------------------------%
    
end
%--------------------------------------------%
% MERGING ALL THE SUBDUCTIONS IN ONE (i.e. a point is PS if it beongs to either subduction zone)


iPS=unique(cell2mat(iPS_sep(:)));
iPS_eff=unique(cell2mat(iPS_eff_sep(:)));

%---- Find BS points ------------------------%
iBS=setdiff(1:Ntetra,iPS_eff);
iBS_moho=find(-pts(:,3)<0 & (depth_moho_pts+pts(:,3))<=0);
iBS_eff=intersect(iBS,iBS_moho);

%--------------------------------------------%


%---- Volumes Statistics --------------------%
i_eff=find(-pts(:,3)<0); % baricenters below zero-level depth
Vol_eff=sum(volume(i_eff,1)); % Ellipsoid Volume effective
Vol_eff_perc=(Vol_eff/Vol_tetra)*100; % Percentage Ellipsoid Volume effective

Vol_PS_eff=sum(volume(iPS_eff,1)); % Ellipsoid Volume PS effective
Vol_BS_eff=sum(volume(iBS_eff,1)); % Ellipsoid Volume BS effective
Vol_PS_BS_eff=Vol_PS_eff+Vol_BS_eff; % Ellipsoid Volume PS+BS effective (smaller than Vol_eff)
Vol_PS_eff_perc=(Vol_PS_eff/Vol_PS_BS_eff)*100;
Vol_BS_eff_perc=(Vol_BS_eff/Vol_PS_BS_eff)*100;
try
   i_PS_BS_eff=[iPS_eff;iBS_eff'];  
catch
   i_PS_BS_eff=[iPS_eff;iBS_eff];   % Removed ' (#Jacopo 11/12/2019)
end
%--------------------------------------------%


%---- PS/BS Lambda --------------------------%
gauss_PS_BS_eff=mvnpdf([ptsx(i_PS_BS_eff),ptsy(i_PS_BS_eff),pts(i_PS_BS_eff,3)*1000],...
    [hypox,hypoy,hypoz*1000],cov); % 3D Gaussian ditsribution for effective PS and BS points
gauss_PS_eff=mvnpdf([ptsx(iPS_eff),ptsy(iPS_eff),pts(iPS_eff,3)*1000],...
    [hypox,hypoy,hypoz*1000],cov); % 3D Gaussian ditsribution for effective PS points
gauss_BS_eff=mvnpdf([ptsx(iBS_eff),ptsy(iBS_eff),pts(iBS_eff,3)*1000],...
    [hypox,hypoy,hypoz*1000],cov); % 3D Gaussian ditsribution for effective BS points

lambda_PS=(sum(gauss_PS_eff.*volume(iPS_eff,1)))/(sum(gauss_PS_BS_eff.*volume(i_PS_BS_eff,1))); % lambda PS
lambda_BS=(sum(gauss_BS_eff.*volume(iBS_eff,1)))/(sum(gauss_PS_BS_eff.*volume(i_PS_BS_eff,1))); % lambda BS

%---- PS Lambda for each zone ---------------% (added #Fabrizio 19/12/2019)
gauss_PS_eff_zone=cell(Nz,1);
lambda_PS_zone=zeros(Nz,1);
for i=1:Nz
    gauss_PS_eff_zone{i,1}=mvnpdf([ptsx(iPS_eff_sep{i,1}),ptsy(iPS_eff_sep{i,1}),pts(iPS_eff_sep{i,1},3)*1000],...
        [hypox,hypoy,hypoz*1000],cov); % 3D Gaussian ditsribution for effective PS points for a single zone
    
    iPS_no_unique=cell2mat(iPS_eff_sep);
    
    gauss_PS_no_unique=mvnpdf([ptsx(iPS_no_unique),ptsy(iPS_no_unique),pts(iPS_no_unique,3)*1000],...
        [hypox,hypoy,hypoz*1000],cov); % 3D Gaussian ditsribution for PS points w/o considering the uniqueness
    
    lambda_PS_zone(i,1)=(sum(gauss_PS_eff_zone{i,1}.*volume(iPS_eff_sep{i,1},1)))/(sum(gauss_PS_no_unique.*volume(iPS_no_unique,1)))*lambda_PS; % lambda PS for each zone

    if (isnan(lambda_PS_zone(i,1)))
        lambda_PS_zone(i,1)=0;
    end
end
    


hypo_deep=0; % flag to plot (1) or not (2) Moho if hypocenter is too deep
if (isnan(lambda_BS) | isnan(lambda_PS))
    hypo_deep=1; % flag to plot Moho if hypocenter is too deep
    
    dist_pts_moho=zeros(Ntetra,1);
    for i=1:Ntetra
        dist_pts_moho(i,1)=min(sqrt((ptsx(i,1)-moho_pts(2).moho_par(:,1)).^2.+(ptsy(i,1)-moho_pts(3).moho_par(:,1)).^2.+(pts(i,3)*1000+moho_pts(4).moho_par(:,1)*1000).^2)/1000);
    end
    [min_pts_moho, imin_pts_moho]=min(dist_pts_moho);
    
    [min_pts_Sub, imin_pts_Sub]=min(cell2mat(distances_sep));
    
    if (min_pts_moho < min_pts_Sub-tsh)
        lambda_BS=1;
        lambda_PS=0;
    else
        lambda_BS=0;
        lambda_PS=1;
    end
end
%--------------------------------------------%


%% ---- Plot Subduction and Ellipsoid ---------%
if false
    figure
    hold on

    for i=1:Nz

        mesh2.faces=mesh_subd(4).mesh_par{i,1}(:,2:end);

        mesh2.vertices=...
            [mesh_subd(3).mesh_par{i,1}(:,2),mesh_subd(3).mesh_par{i,1}(:,3),...
            mesh_subd(3).mesh_par{i,1}(:,4)/1000];

        mesh2_up=mesh2;
        mesh2_up.vertices=mesh2.vertices;
        mesh2_up.vertices(:,3)=mesh2_up.vertices(:,3)+tsh;

        mesh2_down=mesh2;
        mesh2_down.vertices=mesh2.vertices;
        mesh2_down.vertices(:,3)=mesh2_down.vertices(:,3)-tsh;

        [xplane,yplane]=meshgrid(min(pts(:,1))-2:1:max(pts(:,1))+2,min(pts(:,2))-2:1:max(pts(:,2))+2);
        zplane=zeros(size(xplane));


        az=36; % azimuth for 3D view
        el=3; % elevation for 3D view

        subplot(1,2,1)
        patch(mesh2,'FaceAlpha',.3,'EdgeColor','none'); xlabel('Lon'); ylabel('Lat'); zlabel('Depth (km)'); hold on
        patch(mesh2_up,'FaceAlpha',.1,'EdgeColor','none'); xlabel('Lon'); ylabel('Lat'); zlabel('Depth (km)'); hold on
        patch(mesh2_down,'FaceAlpha',.1,'EdgeColor','none'); xlabel('Lon'); ylabel('Lat'); zlabel('Depth (km)'); hold on

        subplot(1,2,2)
        patch(mesh2,'FaceAlpha',.3,'EdgeColor','none'); xlabel('Lon'); ylabel('Lat'); zlabel('Depth (km)'); hold on
        patch(mesh2_up,'FaceAlpha',.1,'EdgeColor','none'); xlabel('Lon'); ylabel('Lat'); zlabel('Depth (km)'); hold on
        patch(mesh2_down,'FaceAlpha',.1,'EdgeColor','none'); xlabel('Lon'); ylabel('Lat'); zlabel('Depth (km)'); hold on


    end

    subplot(1,2,1)
    surf(xplane,yplane,zplane)


    subplot(1,2,2)

    surf(xplane,yplane,zplane)

    if (hypo_deep==1) % i.e. hypocenter too deep
        plot3(pts(iPS,1),pts(iPS,2),-pts(iPS,3),'*r')
        plot3(pts(iBS,1),pts(iBS,2),-pts(iBS,3),'*b')
        plot3(moho_pts(1).moho_par(:,1),moho_pts(1).moho_par(:,2),moho_pts(1).moho_par(:,3),'ok','MarkerSize',1)
        if (lambda_BS==1)
            plot3(pts(imin_pts_moho,1),pts(imin_pts_moho,2),-pts(imin_pts_moho,3),'ob','MarkerFaceColor','b','MarkerEdgeColor','cyan','MarkerSize',10)
            title('Ellipsoid too deep: BS more likely')
        elseif (lambda_PS==1)
            plot3(pts(imin_pts_Sub,1),pts(imin_pts_Sub,2),-pts(imin_pts_Sub,3),'or','MarkerFaceColor','r','MarkerEdgeColor','magenta','MarkerSize',10)
            title('Ellipsoid too deep: PS more likely')
        end
    else
        plot3(pts(iPS_eff,1),pts(iPS_eff,2),-pts(iPS_eff,3),'*r')
        plot3(pts(iBS_eff,1),pts(iBS_eff,2),-pts(iBS_eff,3),'*b')
        title('Ellipsoid, effective PS(red), BS(blue)')
    end
    view(az,el)
    grid on
    axis square
end
%--------------------------------------------%

%======================================================================%
end

