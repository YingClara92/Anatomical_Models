function [i_label,j_label,length_point] = cal_wepl_z(info,roi_center,gantry_angle,couch_angle,Length,ct_orig,RSP,energy)

Nx=size(ct_orig,1)+1; % number of planes in x-axis
Ny=size(ct_orig,2)+1;
Nz=size(ct_orig,3)+1;
dx=info.raw.pixdim(2); % length between each planes in x axis (mm)
dy=info.raw.pixdim(3);
dz=info.raw.pixdim(4);
total_length = Length;
if gantry_angle > 180
    gantry_angle = gantry_angle - 360;
end

spot_coord = []; 
%% 5) WEPL
for spot_indice = 1:(length(spot_coord)+1)
    
    if spot_indice == length(spot_coord)+1 % use the roi center as the reference
%         xpos=Nx-1-roi_center(1);
        xpos = roi_center(1);
        ypos = roi_center(2);
        zpos = roi_center(3);
    else
        xpos=Nx-1-spot_coord(spot_indice,1); % set relative position of the pixel region and the orgin
        ypos=spot_coord(spot_indice,2);
        zpos=spot_coord(spot_indice,3);
    end
    
    bx=-xpos*dx;
    by=-ypos*dy;
    bz=-zpos*dz;
    
    %-----------------------
    
    beam_angle_test=0; % to differentiate whether the beam angle pass through the chest
    
    % convert the linac angle to spherical coordinate angle
    %----------------------
    
    if gantry_angle==180 || gantry_angle==90 || gantry_angle==0 % avoid bug
        gantry_angle = gantry_angle-0.0001;
    end
    if gantry_angle==-180 || gantry_angle==-90 % avoid bug
        gantry_angle = gantry_angle+0.0001;
    end
    
    if couch_angle==90 || couch_angle==0 % avoid bug
        couch_angle = couch_angle-0.0001;
    end
    if couch_angle==-90 % avoid bug
        couch_angle = couch_angle+0.0001;
    end
    
    aa0 = pi*gantry_angle/180;
    alpha = pi*couch_angle/180;
    
    pa = acos(sin(aa0)*sin(alpha))*180/pi; % convert the ga and ca to polar angle
    aa = atan(tan(aa0)*cos(alpha))*180/pi; % convert the ga and ca to azimuthal angle
    
    if abs(gantry_angle) > 90
        aa = aa+180;
    end
    
    test_angle_pa=pa;
    test_angle_aa=180-aa;
    
    pa=pi*((test_angle_pa)/180); % convert the angle from degree to radian
    aa=pi*((test_angle_aa)/180);
    %----------------------
    
    
    
    p1x=(total_length)*sin(pa)*cos(aa); % compute the x coordinate of the starting point
    
    p1y=(total_length)*sin(pa)*sin(aa);
    
    p1z=(total_length)*cos(pa);
    
    p2x=0; % assume the end point of the beam is on the orgin of the coordinate
    p2y=0;
    p2z=0;
    
    
    %find a_min and a_max
    %----------------------
    ax=((bx+dx*(0:(Nx-1)))-p1x)/(p2x-p1x); % compute the a-parameter for each x planes
    ay=((by+dy*(0:(Ny-1)))-p1y)/(p2y-p1y);
    az=((bz+dz*(0:(Nz-1)))-p1z)/(p2z-p1z);
    
    ax_min=min(ax(1),ax(Nx)); % compute the min value of ax
    ax_max=max(ax(1),ax(Nx));
    ay_min=min(ay(1),ay(Ny));
    ay_max=max(ay(1),ay(Ny));
    az_min=min(az(1),az(Nz));
    az_max=max(az(1),az(Nz));
    
    a_min=max(ax_min,max(ay_min,az_min)); % compute a_min
    a_max=min(1,min(ax_max,min(ay_max,az_max)));
    
    pa_minx=p1x+a_min*(p2x-p1x); % compute the x coordinate of a_min
    pa_miny=p1y+a_min*(p2y-p1y);
    pa_minz=p1z+a_min*(p2z-p1z);
    pa_maxx=p1x+a_max*(p2x-p1x);
    pa_maxy=p1y+a_max*(p2y-p1y);
    pa_maxz=p1z+a_max*(p2z-p1z);
    %----------------------
    
    
    % determine i_min, i_max, j_min, j_max
    %----------------------
    % i_min=min(i_f,I_l) and i_max=max((i_f,I_l),
    %where i_f is the number of the first intersected x-plane and i_l is that of the last intersected x-plane
    
    if p1x<p2x  % the beam come from -x to x
        
        if  abs(a_min-ax_min)<=2*eps(a_min) % if a_min=ax_min [There are problem when a_min=ax_min=ay_min]
            i_min=1;
        else
            i_min=ceil((pa_minx-bx)/dx); % [There are problem when i_f & i_l are between the same pair of x planes]
        end
        
        if  abs(a_max-ax_max)<=2*eps(a_max) % if a_max=ax_max
            i_max=Nx-1;
        else
            i_max=floor((pa_maxx-bx)/dx);
        end
        
    end
    
    if p1x>p2x  % the beam come from x to -x
        
        if  abs(a_min-ax_min)<=2*eps(a_min)
            i_max=Nx-2;
        else
            i_max=floor((pa_minx-bx)/dx);
        end
        
        if  abs(a_max-ax_max)<=2*eps(a_max)
            i_min=0;
        else
            i_min=ceil((pa_maxx-bx)/dx);
        end
        
    end
    
    if p1y<p2y  % the beam come from -y to y
        
        if  abs(a_min-ay_min)<=2*eps(a_min)
            j_min=1;
        else
            j_min=ceil((pa_miny-by)/dy);
        end
        
        if  abs(a_max-ay_max)<=2*eps(a_max)
            j_max=Ny-1;
        else
            j_max=floor((pa_maxy-by)/dy);
        end
        
    end
    
    if p1y>p2y  % the beam come from y to -y
        
        if  abs(a_min-ay_min)<=2*eps(a_min)
            j_max=Ny-2;
        else
            j_max=floor((pa_miny-by)/dy);
        end
        
        if  abs(a_max-ay_max)<=2*eps(a_max)
            j_min=0;
        else
            j_min=ceil((pa_maxy-by)/dy);
        end
        
    end
    
    if p1z<p2z  % the beam come from -z to z
        
        if  abs(a_min-az_min)<=2*eps(a_min)
            k_min=1;
        else
            k_min=ceil((pa_minz-bz)/dz);
        end
        
        if  abs(a_max-az_max)<=2*eps(a_max)
            k_max=Nz-1;
        else
            k_max=floor((pa_maxz-bz)/dz);
        end
        
    end
    
    if p1z>p2z  % the beam come from z to -z
        
        if  abs(a_min-az_min)<=2*eps(a_min)
            k_max=Nz-2;
        else
            k_max=floor((pa_minz-bz)/dz);
        end
        
        if  abs(a_max-az_max)<=2*eps(a_max)
            k_min=0;
        else
            k_min=ceil((pa_maxz-bz)/dz);
        end
        
    end
    %----------------------
    
    
    % compute the path length in each pixels
    %----------------------
    ax_i_min=((bx+dx*i_min)-p1x)/(p2x-p1x); % compute the a-parameter of i_min
    ax_i_max=((bx+dx*i_max)-p1x)/(p2x-p1x);
    ay_j_min=((by+dy*j_min)-p1y)/(p2y-p1y);
    ay_j_max=((by+dy*j_max)-p1y)/(p2y-p1y);
    az_k_min=((bz+dz*k_min)-p1z)/(p2z-p1z);
    az_k_max=((bz+dz*k_max)-p1z)/(p2z-p1z);
    
    % determine the a-parameter of the first intersection point of the beam
    % with x,y,z-plane, after the beam enter the voxel space
    
    
    ax0=min(ax_i_min,ax_i_max);
    ay0=min(ay_j_min,ay_j_max);
    az0=min(az_k_min,az_k_max);
    
    % for special case
    %----------------------
    if dx>(abs(pa_minx-p2x))
        ax0=10;
    end
    
    if dy>(abs(pa_miny-p2y))
        ay0=10;
    end
    
    if dz>(abs(pa_minz-p2z))
        az0=10;
    end
    %----------------------
    
    
    % determine the indices of the first voxel after the beam enter the pixel space
    i=floor(((p1x+((min(ax0,min(ay0,az0))+a_min)/2)*(p2x-p1x))-bx)/dx);
    j=floor(((p1y+((min(ax0,min(ay0,az0))+a_min)/2)*(p2y-p1y))-by)/dy);
    k=floor(((p1z+((min(ax0,min(ay0,az0))+a_min)/2)*(p2z-p1z))-bz)/dz);
    
    
    % set the initial parameters for the loop
    Np=(i_max-i_min+1)+(j_max-j_min+1)+(k_max-k_min+1); % total calculating times for the loop
    ac=a_min; % set the intersection point at the back as a_min initially
    axf=ax0; % set the intersection point at the front for x axis as ax0 initally
    ayf=ay0;
    azf=az0;
    plength=zeros(Nx-1,Ny-1,Nz-1); % create a matrix to store the pathlength of each voxels
    
    
    % for special case
    %----------------------
    if dx>(abs(pa_minx-p2x))
        
        Np=(j_max-j_min+1)+(k_max-k_min+1);
        axf=100;
        
        i=i_min-1;
        
    end
    
    if dy>(abs(pa_miny-p2y))
        
        Np=(i_max-i_min+1)+(k_max-k_min+1);
        ayf=100;
        
        j=j_min-1;
        
    end
    
    if dz>(abs(pa_minz-p2z))
        
        Np=(i_max-i_min+1)+(j_max-j_min+1);
        azf=100;
        
        k=k_min-1;
        
    end
    %----------------------
    
    
    % compute the parameters which will be used in the loop
    %----------------------
    axu=dx/abs(p2x-p1x); % compute the increment of ax
    ayu=dy/abs(p2y-p1y);
    azu=dz/abs(p2z-p1z);
    
    if p1x<p2x % determine the increment of x indice
        iu=1;
    else
        iu=-1;
    end
    
    if p1y<p2y
        ju=1;
    else
        ju=-1;
    end
    
    if p1z<p2z
        ku=1;
    else
        ku=-1;
    end
    
    tWEPL=0;
    tlength=0;
    WEPL_3D=zeros(Nx-1,Ny-1,Nz-1);
    soft_WEPL=0;
    bone_WEPL=0;
    
    i0_start = i;
    j0_start = j;
    k0_start = k;
    
    if ct_orig(i+1,j+1,k+1)>-700 % exclude the beam's path which start at patient's body
        beam_angle_test=1;
    end
    flag = 1;
    %----------------------
    
    % loop for computing WEPL
    %----------------------
    i_label = [];
    j_label = [];
    n = 0;
    while (flag == 1)
        n = n+1;
        if i<0 || j<0 || k<0 || i>(Nx-2) || j>(Ny-2) || k>(Nz-2) % stop the loop if the beam exited the voxel space
            break
        end
        if tWEPL > 180
            see = 1;
        end
        if ((axf<ayf) && (axf<azf)) || ((axf<azf) && (axf==ayf))
            plength(i+1,j+1,k+1)=(axf-ac)*total_length; % plength is the path length of that pixel
            tWEPL=tWEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1); % tWEPL is the total WEPL of that point
            tlength=tlength+plength(i+1,j+1,k+1); % tlength is the total path length of that point
            if  tWEPL > energy - 0.1
                i_label = i+1;
                j_label = j+1;
                break;
            end
            if ct_orig(i+1,j+1,k+1)>=900 % pixel is regarded as bone if HU>=900
                bone_WEPL=bone_WEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            else
                soft_WEPL=soft_WEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            end
            i=i+iu;
            ac=axf;
            axf=axf+axu;
            
        elseif (ayf<axf) && (ayf<azf)
            plength(i+1,j+1,k+1)=(ayf-ac)*total_length;
            tWEPL=tWEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            tlength=tlength+plength(i+1,j+1,k+1);
            if  tWEPL > energy - 0.1
                i_label = i+1;
                j_label = j+1;
                break;
            end
            if ct_orig(i+1,j+1,k+1)>=900
                bone_WEPL=bone_WEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            else
                soft_WEPL=soft_WEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            end
            j=j+ju;
            ac=ayf;
            ayf=ayf+ayu;
            
        else
            plength(i+1,j+1,k+1)=(azf-ac)*total_length;
            tWEPL=tWEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            tlength=tlength+plength(i+1,j+1,k+1);
            if  tWEPL > energy - 0.1
                i_label = i+1;
                j_label = j+1;
                break;
            end
            if ct_orig(i+1,j+1,k+1)>=900
                bone_WEPL=bone_WEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            else
                soft_WEPL=soft_WEPL+plength(i+1,j+1,k+1)*RSP(i+1,j+1,k+1);
            end
            k=k+ku;
            ac=azf;
            azf=azf+azu;
        end
        
    end 
end

if isempty(i_label) || isempty(j_label)
    length_point = 0;
else
%     length_point = sqrt(((i+1 - i0_start)*dx)^2+((j+1 - j0_start)*dy)^2);
    length_point = sqrt(((i_label - xpos)*dx)^2+((j_label - ypos)*dy)^2);
    if n > Np
        length_point = - length_point;
    end
end
