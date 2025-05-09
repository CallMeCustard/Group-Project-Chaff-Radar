function eff = effectiveness_score(occlusion, rcs_coh, rcs_inc, w1, w2)
% EFFECTIVENESS_SCORE Combines occlusion and inverse RCS into a single score.
%
% eff = effectiveness_score(occlusion, rcs_coh, rcs_inc, w1, w2)
% Inputs:
%   - occlusion: scalar [0–1], geometric LOS blocking
%   - rcs_coh, rcs_inc: scalar, unnormalized coherent/incoherent RCS
%   - w1, w2: weights for occlusion and RCS terms
%
% Output:
%   - eff: scalar effectiveness score [0–1], higher is better

    rcs_combined = 0.5 * (rcs_coh + rcs_inc);
    rcs_score = 1 - rcs_combined;  % Lower RCS → higher score

    eff = w1 * occlusion + w2 * rcs_score;

    % Clamp between 0 and 1
    eff = max(0, min(1, eff));
end
