function results = exp08_satellite_orbits(makePlots,animate)
if nargin < 1, makePlots = true; end
if nargin < 2, animate = true; end
aEarth = 6378.137;
bEarth = aEarth*(1-1/298.257223563);
mu = 398600.4418;
omegaEarth = 7.2921150e-5;
names = {'LEO','MEO','GEO'};
altitude = [500 20200 35786];
a = aEarth+altitude;
ecc = [0.001 0.010 0.0001];
inc = [51.6 55 0]*pi/180;
raan = [0 0 0];
argPeri = [0 0 0];
M0 = [0 45 0]*pi/180;
n = sqrt(mu./a.^3);
period = 2*pi./n;
summary = table(names',altitude',a',period'/3600,ecc', ...
    'VariableNames',{'Orbit','NominalAltitude_km','SemimajorAxis_km','Period_hours','Eccentricity'});
disp('Experiment 8: two-body Keplerian orbits; Earth rotation angle is zero at t = 0'); disp(summary);
time = linspace(0,24*3600,3000);
tracks = cell(1,3);
colors = lines(3);
if makePlots
    orbitFig = figure('Name','Experiment 8: Inertial Orbits','Color','w');
    orbitAx = axes; hold(orbitAx,'on');
    [sx,sy,sz] = sphere(50);
    surf(orbitAx,aEarth*sx,aEarth*sy,bEarth*sz,'FaceColor',[0.15 0.4 0.7],'EdgeColor','none','HandleVisibility','off');
    satellites = gobjects(1,3);
end
for j = 1:3
    r = orbit_position(a(j),ecc(j),inc(j),raan(j),argPeri(j),M0(j)+n(j)*time);
    angle = omegaEarth*time;
    x = cos(angle).*r(1,:) + sin(angle).*r(2,:);
    y = -sin(angle).*r(1,:) + cos(angle).*r(2,:);
    z = r(3,:);
    lon = atan2d(y,x);
    lat = orbit_latitude(x,y,z,aEarth,bEarth);
    radius = sqrt(sum(r.^2,1));
    assert(all(radius>=a(j)*(1-ecc(j))-1e-7) && all(radius<=a(j)*(1+ecc(j))+1e-7));
    assert(all(abs(lat)<=90) && all(abs(lon)<=180));
    tracks{j} = struct('longitude_deg',lon,'latitude_deg',lat,'radius_km',radius);
    if makePlots
        path = orbit_position(a(j),ecc(j),inc(j),raan(j),argPeri(j),linspace(0,2*pi,700));
        plot3(orbitAx,path(1,:),path(2,:),path(3,:),'Color',colors(j,:),'DisplayName',names{j});
        satellites(j) = plot3(orbitAx,r(1,1),r(2,1),r(3,1),'o','MarkerFaceColor',colors(j,:),'Color',colors(j,:),'HandleVisibility','off');
    end
end
if makePlots
    axis(orbitAx,'equal'); grid(orbitAx,'on'); view(orbitAx,35,24);
    xlabel(orbitAx,'Inertial X (km)'); ylabel(orbitAx,'Inertial Y (km)'); zlabel(orbitAx,'Inertial Z (km)');
    title(orbitAx,'Keplerian orbits'); legend(orbitAx,'Location','best');
    figure('Name','Experiment 8: Ground Tracks','Color','w');
    for j = 1:3
        lon = tracks{j}.longitude_deg; lat = tracks{j}.latitude_deg;
        seam = find(abs(diff(lon))>180)+1;
        lon(seam) = NaN; lat(seam) = NaN;
        subplot(3,1,j); plot(lon,lat,'Color',colors(j,:)); grid on;
        xlim([-180 180]); ylim([-90 90]); xlabel('Longitude (degrees)'); ylabel('Geodetic latitude (degrees)');
        title([names{j} ' ground track over 24 hours']);
    end
    if animate
        for t = linspace(0,24*3600,360)
            if ~ishghandle(orbitFig), break; end
            for j = 1:3
                r = orbit_position(a(j),ecc(j),inc(j),raan(j),argPeri(j),M0(j)+n(j)*t);
                set(satellites(j),'XData',r(1),'YData',r(2),'ZData',r(3));
            end
            drawnow; pause(0.005);
        end
    end
end
results = struct('summary',summary,'time_s',time,'groundTracks',{tracks});
end

function r = orbit_position(a,e,i,raan,argPeri,M)
M = mod(M+pi,2*pi)-pi;
E = M;
for iteration = 1:30
    delta = (E-e*sin(E)-M)./(1-e*cos(E));
    E = E-delta;
    if max(abs(delta(:)))<1e-13, break; end
end
assert(max(abs(E-e*sin(E)-M))<1e-11);
rp = [a*(cos(E)-e);a*sqrt(1-e^2)*sin(E);zeros(size(E))];
Rz1 = [cos(raan) -sin(raan) 0;sin(raan) cos(raan) 0;0 0 1];
Rx = [1 0 0;0 cos(i) -sin(i);0 sin(i) cos(i)];
Rz2 = [cos(argPeri) -sin(argPeri) 0;sin(argPeri) cos(argPeri) 0;0 0 1];
r = Rz1*Rx*Rz2*rp;
end

function latitude = orbit_latitude(x,y,z,a,b)
e2 = 1-b^2/a^2;
p = hypot(x,y);
lat = atan2(z,p*(1-e2));
for iteration = 1:12
    N = a./sqrt(1-e2*sin(lat).^2);
    next = atan2(z+e2*N.*sin(lat),p);
    change = max(abs(next(:)-lat(:)));
    lat = next;
    if change<1e-13, break; end
end
latitude = lat*180/pi;
end
