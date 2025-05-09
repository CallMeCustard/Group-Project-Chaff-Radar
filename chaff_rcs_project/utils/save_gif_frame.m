function gif_path = save_gif_frame(im, gifname_base, t, delay)
% Save a frame to a timestamped GIF in /figures, appending each frame
% Call once per frame (inside loop)

    if nargin < 4
        delay = 0.1;
    end

    % Find top-level figures folder
    root = pwd;
    while ~exist(fullfile(root, 'figures'), 'dir')
        root = fileparts(root);
        if isempty(root)
            error('Could not find figures/ folder.');
        end
    end
    figdir = fullfile(root, 'figures');
    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end

    % Use consistent timestamp across frames
    persistent gif_path_static
    if isempty(gif_path_static)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        gif_path_static = fullfile(figdir, [gifname_base '_' timestamp '.gif']);
    end

    [imind, cm] = rgb2ind(im, 256);

    if t == 1
        imwrite(imind, cm, gif_path_static, 'gif', 'Loopcount', inf, 'DelayTime', delay);
    else
        imwrite(imind, cm, gif_path_static, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end

    gif_path = gif_path_static;
    
    fprintf("Saving GIF frame to: %s\n", gif_path_static);

end
