function save_figure_to_figures(filename_base, ext)
% Save current figure to top-level /figures folder with timestamp
% Usage: save_figure_to_figures('occlusion_scene', 'png')

    if nargin < 2
        ext = 'png';  % default to PNG
    end

    % Navigate up until we find project root with 'figures/' folder
    root = pwd;
    while ~exist(fullfile(root, 'figures'), 'dir')
        root = fileparts(root);
        if isempty(root)
            error('Could not find figures/ folder in any parent directory.');
        end
    end

    % Create folder if needed
    figdir = fullfile(root, 'figures');
    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end

    % Timestamped file path
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    full_filename = fullfile(figdir, [filename_base '_' timestamp '.' ext]);

    % Save current figure
    saveas(gcf, full_filename);
    fprintf('✅ Saved figure to: %s\n', full_filename);
end
