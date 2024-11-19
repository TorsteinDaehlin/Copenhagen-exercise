function relative_cop = FindRelativeCop(dynamic_markers, dynamic_lcs, segments, grf, nof)

sides = fieldnames(grf);

for side = 1:length(sides)
    origo = 0.5*(dynamic_markers.([upper(sides{side}(1)) 'MTH1']) + dynamic_markers.([upper(sides{side}(1)) 'MTH5']));
    vec = grf.(sides{side}).cop - origo;

    for frame = 1:nof
        if isequal(lower(sides{side}),'left')
            relative_cop.(sides{side})(frame,1) = (dot(vec(frame,:),[dynamic_lcs.(['foot_' sides{side}(1)]).epx(frame,1:2) 0])/ ...
                segments.(['foot_' sides{side}(1)]).dist_rad) * -100;
        else
            relative_cop.(sides{side})(frame,1) = (dot(vec(frame,:),[dynamic_lcs.(['foot_' sides{side}(1)]).epx(frame,1:2) 0])/ ...
                segments.(['foot_' sides{side}(1)]).dist_rad) * 100;
        end
        relative_cop.(sides{side})(frame,2) = (dot(vec(frame,:),[dynamic_lcs.(['foot_' sides{side}(1)]).epz(frame,1:2) 0])/ ...
            segments.(['foot_' sides{side}(1)]).len) * 100;
    end
end

end