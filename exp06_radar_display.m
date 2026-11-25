function results = exp06_radar_display(makePlots,animate)
if nargin < 1, makePlots = true; end
if nargin < 2, animate = true; end
savedRng = rng;
rngCleanup = onCleanup(@() rng(savedRng));
rng(6,'twister');
ranges = [30 75 50];
azimuth = [45 120 300];
amplitude = [0.8 0.5 0.9];
vr = [-0.05 -0.08 0.04];
vaz = [0.1 -0.15 0.2];
maxRange = 100;
step = 2;
numFrames = 5*360/step;
x = ranges.*cosd(azimuth);
y = ranges.*sind(azimuth);
initial = table((1:3)',ranges',azimuth',x',y', ...
    'VariableNames',{'Target','Range_km','Azimuth_deg','X_km','Y_km'});
disp('Experiment 6: azimuth is counterclockwise from +X; motion is per frame'); disp(initial);
rangeGrid = 0:0.1:maxRange;
detectionCount = zeros(1,3);
timer = zeros(1,3);
if makePlots
    fig = figure('Name','Experiment 6: Radar Displays','Color','w');
    ax1 = subplot(2,2,1); aLine = plot(rangeGrid,zeros(size(rangeGrid))); grid on;
    xlim([0 maxRange]); ylim([-0.2 1.1]); xlabel('Range (km)'); ylabel('Amplitude'); title('A-scan at current beam angle');
    ax2 = subplot(2,2,3); bPoints = scatter(azimuth,ranges,60,amplitude,'filled'); grid on;
    xlim([0 360]); ylim([0 maxRange]); xlabel('Azimuth (degrees)'); ylabel('Range (km)'); title('B-scan: simulated target positions');
    ax3 = subplot(2,2,[2 4]); hold on; axis equal; grid on;
    theta = linspace(0,2*pi,361);
    for radius = 20:20:maxRange
        plot(radius*cos(theta),radius*sin(theta),'Color',[0.75 0.75 0.75]);
    end
    xlim([-maxRange maxRange]); ylim([-maxRange maxRange]); xlabel('X (km)'); ylabel('Y (km)');
    title('PPI: red truth, green recent detections');
    sweepLine = plot([0 maxRange],[0 0],'g-','LineWidth',1.5);
    truthPoints = scatter(x,y,35,'r','filled');
    hitPoints = scatter(nan(1,3),nan(1,3),100,'g');
end
for frame = 1:numFrames
    sweep = mod((frame-1)*step,360);
    ranges = ranges + vr;
    above = ranges > maxRange;
    below = ranges < 0;
    ranges(above) = 2*maxRange-ranges(above);
    ranges(below) = -ranges(below);
    vr(above | below) = -vr(above | below);
    azimuth = mod(azimuth+vaz,360);
    separation = abs(mod(azimuth-sweep+180,360)-180);
    hits = separation <= step;
    detectionCount = detectionCount + hits;
    timer = max(timer-1,0);
    timer(hits) = 25;
    x = ranges.*cosd(azimuth);
    y = ranges.*sind(azimuth);
    ascan = 0.025*randn(size(rangeGrid));
    for j = find(hits)
        ascan = ascan + amplitude(j)*exp(-0.5*((rangeGrid-ranges(j))/0.2).^2);
    end
    if makePlots && (animate || frame==numFrames)
        if ~ishghandle(fig), break; end
        set(aLine,'YData',ascan);
        set(bPoints,'XData',azimuth,'YData',ranges);
        set(sweepLine,'XData',[0 maxRange*cosd(sweep)],'YData',[0 maxRange*sind(sweep)]);
        set(truthPoints,'XData',x,'YData',y);
        hx = x; hy = y;
        hx(timer==0) = NaN; hy(timer==0) = NaN;
        set(hitPoints,'XData',hx,'YData',hy);
        title(ax1,sprintf('A-scan at %.0f degrees',sweep));
        drawnow;
        if animate, pause(0.005); end
    end
end
assert(all(ranges>=0 & ranges<=maxRange));
assert(max(abs(hypot(x,y)-ranges))<1e-10);
results = struct('initial',initial,'finalRange_km',ranges,'finalAzimuth_deg',azimuth,'detectionCount',detectionCount,'framesCompleted',frame);
end
