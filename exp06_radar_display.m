rng(6);

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

disp('Experiment 6: Radar Display');
disp(initial);

rangeGrid = 0:0.1:maxRange;

detectionCount = zeros(1,3);
timer = zeros(1,3);

figure('Name','Experiment 6: Radar Displays','Color','w');

subplot(2,2,1);
aLine = plot(rangeGrid,zeros(size(rangeGrid)));
grid on;
xlim([0 maxRange]);
ylim([-0.2 1.1]);
xlabel('Range (km)');
ylabel('Amplitude');
title('A-scan');

subplot(2,2,3);
bPoints = scatter(azimuth,ranges,60,amplitude,'filled');
grid on;
xlim([0 360]);
ylim([0 maxRange]);
xlabel('Azimuth (degrees)');
ylabel('Range (km)');
title('B-scan');

subplot(2,2,[2 4]);
hold on;
axis equal;
grid on;

theta = linspace(0,2*pi,361);

for radius = 20:20:maxRange
    plot(radius*cos(theta),radius*sin(theta),'Color',[0.75 0.75 0.75]);
end

xlim([-maxRange maxRange]);
ylim([-maxRange maxRange]);
xlabel('X (km)');
ylabel('Y (km)');
title('PPI: red truth, green recent detections');

sweepLine = plot([0 maxRange],[0 0],'g-','LineWidth',1.5);
truthPoints = scatter(x,y,35,'r','filled');
hitPoints = scatter(nan(1,3),nan(1,3),100,'g');

for frame = 1:numFrames

    sweep = mod((frame-1)*step,360);

    ranges = ranges + vr;

    for j = 1:3
        if ranges(j) > maxRange
            ranges(j) = 2*maxRange-ranges(j);
            vr(j) = -vr(j);
        elseif ranges(j) < 0
            ranges(j) = -ranges(j);
            vr(j) = -vr(j);
        end
    end

    azimuth = mod(azimuth+vaz,360);

    separation = abs(mod(azimuth-sweep+180,360)-180);
    hits = separation <= step;

    detectionCount = detectionCount + hits;

    timer = max(timer-1,0);
    timer(hits) = 25;

    x = ranges.*cosd(azimuth);
    y = ranges.*sind(azimuth);

    ascan = 0.025*randn(size(rangeGrid));

    for j = 1:3
        if hits(j)
            ascan = ascan + amplitude(j)* ...
                exp(-0.5*((rangeGrid-ranges(j))/0.2).^2);
        end
    end

    set(aLine,'YData',ascan);
    title(sprintf('A-scan at %.0f degrees',sweep));

    set(bPoints,'XData',azimuth,'YData',ranges);

    set(sweepLine,'XData',[0 maxRange*cosd(sweep)], ...
                  'YData',[0 maxRange*sind(sweep)]);

    set(truthPoints,'XData',x,'YData',y);

    hx = x;
    hy = y;

    hx(timer==0) = NaN;
    hy(timer==0) = NaN;

    set(hitPoints,'XData',hx,'YData',hy);

    drawnow;
    pause(0.005);

end

disp('Final Target Ranges (km):');
disp(ranges);

disp('Final Target Azimuths (degrees):');
disp(azimuth);

disp('Total Detections:');
disp(detectionCount);

disp('Frames Completed:');
disp(frame);

assert(all(ranges>=0 & ranges<=maxRange));
assert(max(abs(hypot(x,y)-ranges))<1e-10);
