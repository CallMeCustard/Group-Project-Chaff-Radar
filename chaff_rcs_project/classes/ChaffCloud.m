classdef ChaffCloud
    properties
        positions       % Nx3 matrix of dipole positions
        orientations    % Nx3 matrix of unit vectors
        lengths         % Nx1 vector (optional, for segment-based models)
        modelType       % string label, e.g. 'sixdof', 'mc-gaussian', etc.
        created         % datetime
    end

    methods
        function obj = ChaffCloud(positions, orientations, lengths, modelType)
            % Constructor with defaults
            if nargin < 4, modelType = 'unspecified'; end
            obj.positions = positions;
            obj.orientations = orientations;
            obj.lengths = lengths;
            obj.modelType = modelType;
            obj.created = datetime('now');
        end

        function plot(obj, scalarField)
            % Plot with optional color field
            if nargin < 2, scalarField = []; end
            figure;
            scatter3(obj.positions(:,1), obj.positions(:,2), obj.positions(:,3), ...
                20, scalarField, 'filled');
            if ~isempty(scalarField), colormap(jet); colorbar; end
            hold on;
            quiver3(obj.positions(:,1), obj.positions(:,2), obj.positions(:,3), ...
                obj.orientations(:,1), obj.orientations(:,2), obj.orientations(:,3), ...
                0.5, 'k');
            axis equal; grid on; view(3);
            title(['Chaff Cloud: ', obj.modelType]);
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end

        function n = count(obj)
            n = size(obj.positions, 1);
        end

        function export(obj, filename)
            % Save object to .mat
            if nargin < 2
                filename = ['chaffcloud_', obj.modelType, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat'];
            end
            chaff = obj; %#ok<NASGU>
            save(filename, 'chaff');
        end
    end
end
