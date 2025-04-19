clear;
clc;
clf();

// Subroutine: Generate coordinates    
function [a1, a2, a3] = lattice_pos(index, Lx, Ly, Lz, lattice_const)
    l = 1;  // Start count of particle number
    for i = -Lx/2:Lx/2                      
        x = i * lattice_const;
        for j = -Ly/2:Ly/2
            y = j * lattice_const;
            for k = -Lz/2:Lz/2
                z = k * lattice_const;
                if (abs(x) < Lx/2) & (abs(y) < Ly/2) & (abs(z) < Lz/2) & (l <= index) then
                    a1 = x;
                    a2 = y;
                    a3 = z;
                    l = l + 1;
                end
            end
        end  
    end  
endfunction

// Initialization
function [x, y, z, xm, ym, zm, vx, vy, vz] = initialization(npart, dt, Lx, Ly, Lz, lattice_const, temp)
    sumvx = 0;
    sumvy = 0;
    sumvz = 0;
    sumv2 = 0;
    for i = 1:npart
        [x(i), y(i), z(i)] = lattice_pos(i, Lx, Ly, Lz, lattice_const);
       
        vx(i) = (rand() - 0.5);
        vy(i) = (rand() - 0.5);
        vz(i) = (rand() - 0.5);

        sumvx = sumvx + vx(i);
        sumvy = sumvy + vy(i);
        sumvz = sumvz + vz(i);
       
        sumv2 = sumv2 + (vx(i)^2 + vy(i)^2 + vz(i)^2);
    end
    sumvx = sumvx/npart;
    sumvy = sumvy/npart;
    sumvz = sumvz/npart;
       
    sumv2 = sumv2/npart;
       
    fs = sqrt(3*temp/sumv2);
   
    for i = 1:npart
        vx(i) = (vx(i) - sumvx) * fs;
        vy(i) = (vy(i) - sumvy) * fs;
        vz(i) = (vz(i) - sumvz) * fs;

        xm(i) = x(i) - vx(i)*dt;
        ym(i) = y(i) - vy(i)*dt;
        zm(i) = z(i) - vz(i)*dt;
    end
endfunction

// Force Calculation
function [fx, fy, fz, en] = force(npart, x, y, z, rc, Lx, Ly, Lz, ecut)
    en = 0;
    fx = zeros(1, npart);
    fy = zeros(1, npart);
    fz = zeros(1, npart);
   
    for i = 1:npart-1
        for j = i+1:npart
            xr = x(i) - x(j);
            xr = xr - Lx*round(xr/Lx);
            yr = y(i) - y(j);
            yr = yr - Ly*round(yr/Ly);
            zr = z(i) - z(j);
            zr = zr - Lz*round(zr/Lz);
 
            r2 = xr^2 + yr^2 + zr^2;
 
            rc2 = rc^2;
           
            if (r2 <= rc2) then
                r2i = 1/r2;
                r6i = r2i^3;
 
                ff = 48*r2i*r6i*(r6i-0.5); //Lennard-Jones force
               
                fx(i) = fx(i) + ff*xr;
                fy(i) = fy(i) + ff*yr;
                fz(i) = fz(i) + ff*zr;
               
                fx(j) = fx(j) - ff*xr;
                fy(j) = fy(j) - ff*yr;
                fz(j) = fz(j) - ff*zr;
                en = en + 4*r6i*(r6i-1) - ecut;
            end
        end
    end
endfunction

// Integrate Equations of Motion
function [temperature, etot, xm, ym, zm, x, y, z, KE, PE] = integrate_eqn(npart, Lx, Ly, Lz, dt, x, y, z, xm, ym, zm, fx, fy, fz, en, temp)
    sumvv2 = 0;
    for i = 1:npart
        xx = 2*x(i) - xm(i) + (dt^2)*fx(i);
        yy = 2*y(i) - ym(i) + (dt^2)*fy(i);
        zz = 2*z(i) - zm(i) + (dt^2)*fz(i);
 
        vxx(i) = (xx - xm(i))/(2*dt);
        vyy(i) = (yy - ym(i))/(2*dt);
        vzz(i) = (zz - zm(i))/(2*dt);
 
        sumvv2 = sumvv2 + (vxx(i)^2 + vyy(i)^2 + vzz(i)^2);
       
        xm(i) = x(i);
        x(i) = xx;
        ym(i) = y(i);
        y(i) = yy;
        zm(i) = z(i);
        z(i) = zz;
        
    end
    temperature = sumvv2/(3*npart);
    KE = sumvv2/(2*npart);
    PE = en/npart;
    etot = KE + PE;
endfunction

// Main Parameters
npart = 100;
Lx = 10;
Ly = 10;
Lz = 10;
lattice_const = 1.1;
rc = 2.5;
ecut = 4*((1/rc)^12 - (1/rc)^6);
dt = 0.001;
temp = 2;
tmax = 100;

fd = mopen("traj_100.xyz", "wt");
energy_file = mopen("energy_log.txt", "wt");
mfprintf(energy_file, "Timestep\tKineticEnergy\tPotentialEnergy\tTotalEnergy\n");


[x, y, z, xm, ym, zm, vx, vy, vz] = initialization(npart, dt, Lx, Ly, Lz, lattice_const, temp);

t = 0;
count = 0;

while (t < tmax)
    [fx, fy, fz, en] = force(npart, x, y, z, rc, Lx, Ly, Lz, ecut);
    [temperature, etot, xm, ym, zm, x, y, z, KE, PE] = integrate_eqn(npart, Lx, Ly, Lz, dt, x, y, z, xm, ym, zm, fx, fy, fz, en, temp);
   
    mfprintf(fd, "%d\n\n", npart);
    for i = 1:npart
        mfprintf(fd, "Ar %f %f %f\n", x(i), y(i), z(i));
    end
    mfprintf(energy_file, "%f\t%f\t%f\t%f\n", t, KE, PE, etot);
   
    t = t + dt;
    count = count + 1;
    disp(count);
end

mfprintf(fd, 'END');
mclose(fd);
mclose(energy_file);
