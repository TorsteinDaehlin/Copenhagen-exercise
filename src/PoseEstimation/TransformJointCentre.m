function global_jc = TransformJointCentre(lcs, local_jc, nof)

% Preallocate
global_jc = zeros(nof,3);

for frame = 1:nof
    % Define rotation matrix
    R = [lcs.epx(frame,:)' lcs.epy(frame,:)' lcs.epz(frame,:)'];

    % Transform joint centre to global frame
    global_jc(frame,:) = (R*local_jc' + lcs.origin(frame,:)')';
end

end
