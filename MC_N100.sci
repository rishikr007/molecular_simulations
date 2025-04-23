function [a1,a2,a3] = lattice_pos(index,Lx,Ly,Lz,lattice_const)
   l=1;                                     
   for i = -Lx/2:Lx/2                       
       x=i*lattice_const; 
       for j = -Ly/2:Ly/2
           y=j*lattice_const; 
           for k = -Lz/2:Lz/2
               z=k*lattice_const;
               if (abs(x)<Lx/2) & (abs(y)<Ly/2) & (abs(z)<Lz/2) & (l<=index) then 
                  a1=x; a2=y; a3=z;
                  l=l+1; 
               end
           end
       end
    end   
    
endfunction


function[x,y,z]=initialize(Lx,Ly,Lz,lattice_const)
    for i=1:npart;
        [x(i),y(i),z(i)]=lattice_pos(i,Lx,Ly,Lz,lattice_const)
        mprintf("%.4f %.4f %.4f\n",x(i), y(i), z(i)) //first thing that gets printed
    end
endfunction


function[en]=potential_calculation(npart,x,y,z,rc,Lx,Ly,Lz,ecut)
    en=0
    for i=1:npart-1
        for j=i+1:npart
            xr=x(i)-x(j)
            xr=xr-Lx*round(xr/Lx)
            yr=y(i)-y(j)
            yr=yr-Ly*round(yr/Ly)
            zr=z(i)-z(j)
            zr=zr-Lz*round(zr/Lz)
            r2=xr^2+yr^2+zr^2
            rc2=rc^2
            if (r2<=rc2)then
                r2i=1/r2
                r6i=r2i^3
                en = en+4*r6i*(r6i-1)-ecut
            end
    end
end
endfunction


function[x,y,z,o,enn, acc_ratio] = displace(npart, Lx,Ly,Lz,lattice_const,del,x,y,z)
    
    total_counter = 0
    accept_counter = 0

    while accept_counter < npart
        
        o=int(rand()*npart) + 1
    
        [en] = potential_calculation(npart, x, y, z, rc, Lx, Ly, Lz, ecut) 
        eno = en
    
        temp_x = x(o) 
        temp_y = y(o)
        temp_z = z(o)
    
        x(o) = x(o) + (rand() - 0.5)*del 
        y(o) = y(o) + (rand() - 0.5)*del
        z(o) = z(o) + (rand() - 0.5)*del

        [en]=potential_calculation(npart,x,y,z,rc,Lx,Ly,Lz,ecut)
        enn = en

        if rand() < exp(-b*(enn - eno)) then
            accept_counter = accept_counter + 1
        else
            x(o) = temp_x
            y(o) = temp_y
            z(o) = temp_z
            
            enn = eno
        end
        total_counter = total_counter + 1
    end
    mprintf('Accepted Counter: %d\t\n', accept_counter)
    mprintf('Total Counter: %d\t\n', total_counter)
    acc_ratio = accept_counter/total_counter
    
endfunction


npart = 100;
Lx=10;
Ly=10;
Lz=10;
rc=2.5;
lattice_const=1.1;
ecut=4*((1/rc)^12 - (1/rc)^6);
temp=2;
del=0.25;
b=1;

count=0

fd=mopen('traj_mc_100.xyz','wt')
energy_mc=mopen('energy_mc_100.txt','wt')
mfprintf(fd, '%d\n\n', npart)
[x,y,z] = initialize(Lx, Ly, Lz, lattice_const)

while count < 10000
    [x, y, z, o, enn, acc_ratio] = displace(npart, Lx, Ly, Lz, lattice_const, del, x, y, z)
    
        if acc_ratio < 0.45 then
        //disp('Acceptance ratio is small. Decreasing del')
        del = del*0.95
    elseif acc_ratio > 0.55 then
        //disp('Acceptance ratio is large. Increasing del')
        del = del*1.05
    end

    mprintf('Cycle: %d\t\tEnergy: %.4f\n', count,enn/npart) 
    mfprintf(energy_mc, '%.4f\t %.4f\t\n',count,enn/npart)
    for i=1:npart
        mfprintf(fd,'Ar %.4f %.4f %.4f\n',x(i),y(i),z(i));
    end
    mfprintf(fd, '%d\n\n', npart);
    count=count+1
    disp(count)
end
mfprintf(fd,'END');
mclose(fd);
mclose(energy_mc)

