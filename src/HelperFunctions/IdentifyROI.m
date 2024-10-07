function  [roi] = IdentifyROI(t, grf, mass)

% Truncate force signal
grf_norm = grf ./ (mass * 9.81);

% The vertical ground reaction force is heavily filtered using a moving
% average smoothing approach before identifying when the participant is in
% contact with the force platform. This region is subsequently refined into
% the regions of interest (ROI)
idx = findchangepts(movmean(grf_norm(:, 3), 200), 'MaxNumChanges', 2);
roi = findchangepts(movmean(grf_norm(idx(1):idx(2), 3), 50), 'MaxNumChanges', 6);

% Reshape ideces so the first column contains the start of a region of
% interest and the second constains the end of a region of interst
roi = reshape(idx(1) + roi, [2, 3])';

% Create data check figure
cmap = parula(9);
cmap = flipud([cmap(1,:); cmap(4, :); cmap(7,:)]);

figure();
c = size(grf_norm, 2);
for i = 1:c
    p (i) = plot(t, grf_norm(:, i), 'Color', cmap(i, :));
    hold on;
    xline(t(roi(i, :)),'k--', {'$t_0$', '$t_{end}$'}, ...
        'Interpreter', 'latex');
end
legend([p(1), p(2), p(3)], {'X', 'Y', 'Z'});
xlabel('$Time (s)$', 'Interpreter', 'latex');
ylabel('$GRF (\frac{N}{m \cdot g})$', 'Interpreter', 'latex');
title('$ROI data check$', 'Interpreter', 'latex');

% Save to output directory

end