function combinedDroneUAV(urgency, distanceKm, needsColdChain)
    %% --- Drone Selection Section ---
    droneNames = {'Wingcopter 198', 'Zipline P2 Zip', 'Redwing Labhawk', ...
                  'Matternet M2', 'Sagar MedCOPTER', 'TechEagle X3'};
    droneRanges = [75, 130, 60, 20, 15, 100];
    droneSpeeds = [150, 120, 100, 80, 60, 120];
    droneColdChain = [true, true, true, true, false, true];
    droneUseType = {'emergency', 'emergency', 'routine', 'routine', 'routine', 'emergency'};

    suitableIdx = [];
    for i = 1:length(droneNames)
        if strcmp(droneUseType{i}, urgency) && ...
           droneRanges(i) >= distanceKm && ...
           (~needsColdChain || droneColdChain(i))
            suitableIdx(end+1) = i;
        end
    end

    if isempty(suitableIdx)
        fprintf('No suitable drone found.\n');
        return;
    end

    [~, maxIdx] = max(droneSpeeds(suitableIdx));
    chosenIdx = suitableIdx(maxIdx);

    fprintf('Selected Drone: %s\n', droneNames{chosenIdx});
    fprintf('Range: %d km\n', droneRanges(chosenIdx));
    fprintf('Speed: %d km/h\n', droneSpeeds(chosenIdx));
    fprintf('Cold Chain Support: %s\n\n', mat2str(droneColdChain(chosenIdx)));

    UAV_speed = droneSpeeds(chosenIdx) / 3.6; % m/s

    %% --- Coordinates ---
    lat = [21.24987208329769, 21.25807802569907];
    lon = [81.60473512913642, 81.57956938003011];

    %% --- Manual Elevation Inputs ---
    altStartTerrain = 290;
    altEndTerrain   = 283;
    buildingHeightOffset = 15;

    startIsBuilding = true;
    endIsBuilding = true;

    altStart = altStartTerrain + buildingHeightOffset * startIsBuilding;
    altEnd   = altEndTerrain + buildingHeightOffset * endIsBuilding;

    cruiseAltitudeAGL = 80;
    altCruise = altStart + cruiseAltitudeAGL;

    scale = 1e5;
    start_scaled = [lat(1), lon(1)] * scale;
    end_scaled   = [lat(2), lon(2)] * scale;

    flatDist = norm(end_scaled - start_scaled);
    totalDist = sqrt(flatDist^2 + cruiseAltitudeAGL^2);
    estTime = max(15, totalDist / UAV_speed * 1.5);  % Smoothed time to reduce jerk
    t = linspace(0, estTime, 100);

    % --- Generate trajectory and derivatives ---
    [qx, qx_dot, qx_ddot, qx_jerk] = septicTrajectory(t, estTime, start_scaled(1), end_scaled(1));
    [qy, qy_dot, qy_ddot, qy_jerk] = septicTrajectory(t, estTime, start_scaled(2), end_scaled(2));
    [qz, qz_dot, qz_ddot, qz_jerk] = septicTrajectory(t, estTime, altStart, altCruise);

    qz(qz < altStart) = altStart;

    bestPath = [qx', qy', qz'];
    bestLat = bestPath(:,1) / scale;
    bestLon = bestPath(:,2) / scale;
    bestAltitudes = bestPath(:,3);

    segment_lengths = vecnorm(diff(bestPath), 2, 2);
    totalDist = sum(segment_lengths);
    bestTime = totalDist / UAV_speed;

    %% --- Plot: 3D Globe Path ---
    uif = uifigure;
    g = geoglobe(uif);
    geoplot3(g, lat, lon, [altStart altEnd], 'ro', ...
             'MarkerSize', 10, 'LineWidth', 2, 'HeightReference', 'terrain');
    geoplot3(g, bestLat, bestLon, bestAltitudes, 'b', ...
             'LineWidth', 2.5, 'HeightReference', 'terrain');

    %% --- Plot: Position, Velocity, Acceleration, Jerk ---
    figure('Name', 'Trajectory Derivatives');

    subplot(4,1,1); plot(t, qz, 'b', 'LineWidth', 2); ylabel('Position Z (m)'); title('Z-Direction Motion Profile'); grid on;
    subplot(4,1,2); plot(t, qz_dot, 'r', 'LineWidth', 2); ylabel('Velocity Z (m/s)'); grid on;
    subplot(4,1,3); plot(t, qz_ddot, 'k', 'LineWidth', 2); ylabel('Acceleration Z (m/s^2)'); grid on;
    subplot(4,1,4); plot(t, qz_jerk, 'm', 'LineWidth', 2); ylabel('Jerk Z (m/s^3)'); xlabel('Time (s)'); grid on;

    [~, jerkPeaks] = findpeaks(abs(qz_jerk));
    fprintf("Number of jerk events in Z-direction: %d\n", numel(jerkPeaks));

    %% --- Output Info ---
    fprintf("\nSeptic Trajectory Path\n");
    fprintf("From: [%.6f, %.6f] | Elevation: %.2f m | Final Start: %.2f m\n", ...
        lat(1), lon(1), altStartTerrain, altStart);
    fprintf("To  : [%.6f, %.6f] | Elevation: %.2f m | Final End: %.2f m\n", ...
        lat(2), lon(2), altEndTerrain, altEnd);
    fprintf("Cruise Altitude AGL: %.2f m\n", cruiseAltitudeAGL);
    fprintf("Total 3D Distance: %.2f m\n", totalDist);
    fprintf("Estimated Flight Time: %.2f s\n", bestTime);
end

function [pos, vel, acc, jerk] = septicTrajectory(t, T, p0, pf)
    v0 = 0; vf = 0; a0 = 0; af = 0; j0 = 0; jf = 0;
    A = [1 0 0 0 0 0    0     0;
         0 1 0 0 0 0    0     0;
         0 0 2 0 0 0    0     0;
         0 0 0 6 0 0    0     0;
         1 T T^2 T^3 T^4 T^5 T^6 T^7;
         0 1 2*T 3*T^2 4*T^3 5*T^4 6*T^5 7*T^6;
         0 0 2 6*T 12*T^2 20*T^3 30*T^4 42*T^5;
         0 0 0 6 24*T 60*T^2 120*T^3 210*T^4];

    b = [p0; v0; a0; j0; pf; vf; af; jf];
    coeffs = A \ b;
    a = flip(coeffs');

    pos = polyval(a, t);
    vel = polyval(polyder(a), t);
    acc = polyval(polyder(polyder(a)), t);
    jerk = polyval(polyder(polyder(polyder(a))), t);
end
