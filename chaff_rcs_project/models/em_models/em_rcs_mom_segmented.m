function [rcs_array, total_rcs] = em_rcs_mom_segmented(positions, orientations, lengths, lambda, num_segments)
% EM_RCS_MOM_SEGMENTED
% Simplified Zhou-style MoM RCS estimation using segmented dipoles
% Based on: Zhou et al., IEEE TAP 2023

N = size(positions, 1);
k = 2*pi / lambda;
E_wave = [1 0 0];  % Radar wave from +X

% Prepare outputs
rcs_array = zeros(N,1);

% Loop over dipoles
for i = 1:N
    L = lengths(i);
    ori = orientations(i,:) / norm(orientations(i,:));
    P = positions(i,:);

    seg_len = L / num_segments;
    E_total = 0;

    for s = 1:num_segments
        % Segment center point along dipole
        frac = (s - 0.5) / num_segments - 0.5;
        seg_pos = P + frac * L * ori;

        % Incident phase at this segment
        r_hat = E_wave / norm(E_wave);
        phase = -k * dot(seg_pos, r_hat);

        % Segment current (sinusoidal approx)
        I_seg = sin(pi * abs(frac)); % peak at center, zero at ends
        E_seg = I_seg * seg_len * sin(acos(dot(ori, r_hat))) * exp(1j * phase);

        E_total = E_total + E_seg;
    end

    % RCS ∝ |E|^2 normalized to lambda^2
    rcs_array(i) = (abs(E_total)^2) / lambda^2;
end

total_rcs = sum(rcs_array);
end
