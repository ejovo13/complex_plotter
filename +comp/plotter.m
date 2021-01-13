classdef plotter
    %PLOTTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f              % f(z)
        bounds
    end
    
    properties(Hidden)
        C               % The color values used for contour projections and 2d color plots
        Z               % The matrix of input values for f(z)
        foZ             % The matrix of f(z)
        plotType
    end
    
    methods        
        function obj = plotter(anonFun, bounds, plotType)
            %PLOTTER Construct an instance of this class
            %   Pass an anonymous complex function to the constructor that
            %   you wish to be plotted
            assert(isa(anonFun, 'function_handle'), "First parameter must be a function handle")
            if nargin == 1
                bounds = [-10 10 -10 10];
                plotType = '3d';
            elseif nargin == 2
                plotType = '3d';
            end
            
            obj.f = anonFun;
            obj.plotType = plotType;
            obj.bounds = bounds; % Set the bounds and compute f(z)
            obj.plot3;
        end
        
        %% Plotting methods
        
        function S = plotReal3d(obj) 
            %PLOTREAL3D Plot the real component of f(z)
            %   Using the real and imaginary parts of z as the x and y axis,
            %   plot the real component of f(z) on the z-axis
            S = surfc(real(obj.Z), imag(obj.Z), real(obj.foZ));
            S(2).LevelStep = S(2).LevelStep/10; % Add ten times the number of contours
            xlabel('Re$(z)$', 'Interpreter', 'latex');
            ylabel('Im$(z)$', 'Interpreter', 'latex');
            zlabel('Re$(f(z))$', 'Interpreter', 'latex');
            title("Real part of $f(z)$", 'Interpreter', 'latex');
            shading interp
        end
        
        function plotImag3d(obj)
            %PLOTIMAG3D Plot the imag component of f(z)
            %   Using the real and imaginary parts of z as the x and y axis,
            %   plot the imag component of f(z) on the z-axis
            S = surfc(real(obj.Z), imag(obj.Z), imag(obj.foZ));
            S(2).LevelStep = S(2).LevelStep/10; % Add ten times the number of contours
            xlabel('Re$(z)$', 'Interpreter', 'latex');
            ylabel('Im$(z)$', 'Interpreter', 'latex');
            zlabel('Im$(f(z))$', 'Interpreter', 'latex');
            title("Imaginary part $f(z)$", 'Interpreter', 'latex');
            shading interp
        end
        
        function plotMod3d(obj)
            %PLOTMOD3D Plot the modulus of f(z)
            %   Using the real and imaginary parts of z as the x and y axis,
            %   plot the modulus of f(z) on the z-axis
            S = surfc(real(obj.Z), imag(obj.Z), abs(obj.foZ));
            S(2).LevelStep = S(2).LevelStep/10; % Add ten times the number of contours
            xlabel('Re$(z)$', 'Interpreter', 'latex');
            ylabel('Im$(z)$', 'Interpreter', 'latex');
            zlabel('$|f(z)|$', 'Interpreter', 'latex');
            title("Modulus of $f(z)$", 'Interpreter', 'latex');
            shading interp
        end
        
        function plot3(obj, funHandle)

            if nargin == 2
                obj.f = funHandle;
            end

            %PLOT3 Plot the Real, Imag, and Modulus of f(z) in the same figure
            fun_name = string(char(obj.f));
            fun_name = erase(fun_name, "@(z)");
            figure('NumberTitle', 'off', 'Name', join(["3d plots of" fun_name]));
            % disp(fun_name);
            tiledlayout(1, 4);
            nexttile
            obj.plotReal3d;
            nexttile
            obj.plotImag3d;
            nexttile
            obj.plotMod3d;            
            colormap jet            
            shading interp
            a = nexttile;
            obj.plotColorSimple;
            colormap(a, 'hsv');
            axis(a, 'square')

        end
        
        %% Domain Plotting

        function plotColorSimple(obj)
            dim = size(obj.foZ);
            c_dim = [dim, 3];

            loMod = @(m) 1 - 0.5^m;

            H = angle(obj.foZ);
            S = ones(dim);
            V = arrayfun(loMod, abs(obj.foZ));
            
            Colors = zeros(c_dim);
            Colors(:,:,1) = H;
            Colors(:,:,2) = S;
            Colors(:,:,3) = V;

            P = pcolor(real(obj.Z), imag(obj.Z), H);
            % shading interp
            set(P, 'EdgeColor', 'none');
            xlabel('Re$(z)$', 'Interpreter', 'latex');
            ylabel('Im$(z)$', 'Interpreter', 'latex');
            title("Domain coloring plot of $f(z)$", 'Interpreter', 'latex');
            colorbar
        end

        function plotReimanSurface(obj)

            surf(real(obj.Z), imag(obj.Z), real(obj.foZ), imag(obj.foZ));
            shading interp
            colormap hsv
        end


        function zoom(obj, xy_value) 

            obj.bounds = [-xy_value xy_value -xy_value xy_value];
            obj.plot3;

        end










        function obj = makeGrid(obj, size)
            %MAKEGRID Generate Z, the array of inputs to the function f(z)
                % Generate a 1000x1000 meshed grid of complex values.
                %
                %BOUNDS - For z = x + iy, the bounds shall be passed as
                %         [xmin xmax ymin ymax]
            if nargin < 3
                size = 500;
                % disp('using size 100')
            end
            
                
            X = linspace(obj.bounds(1), obj.bounds(2), size);
            Y = linspace(obj.bounds(3), obj.bounds(4), size);
            
            [X, Y] = meshgrid(X, Y);
            Z_ = X + Y*i;
            
            obj.Z = Z_;         
        end
        
        
        %% Getter and Setter Methods
        
        function obj = setColor
            
        end

        function obj = set.f(obj, newFun) 

            obj.f = newFun;
            obj = obj.applyFun;
            % disp("set.f called");
        end

        function obj = set.bounds(obj, bounds)
            obj.bounds = bounds;
            obj = obj.makeGrid;
            obj = obj.applyFun;
        end

        
        
        %% Helper Methods
        
        function obj = applyFun(obj)         
            obj.foZ = arrayfun(obj.f, obj.Z);            
        end
        
        
        
        
        
        
        
        function obj = plot_input(obj)
            X = obj.grid_lower:obj.step_size:obj.grid_upper;
            Y = X;
            [X, Y] = meshgrid(X, Y);
            Z = X + Y*i;
            obj.C = sqrt(X.^2 + Y.^2);
            pcolor(real(Z), imag(Z), obj.C);            
            obj.Z_in = Z;
            axis equal
        end
        
        function obj = plot_output(obj)
            Z = arrayfun(obj.comp_fun, obj.Z_in);
            pcolor(real(Z), imag(Z), obj.C);
            axis equal
            obj.Z_out = Z;
        end
        
        function plot(obj)
            tiledlayout(1,2);
            nexttile;
            obj = obj.plot_input;
            hold on
            obj.draw_box;
            hold off
            
            ax = gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            
            nexttile
            obj = obj.plot_output;
            hold on
            obj.draw_box;
            hold off
            
            ax = gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            
            axis equal
            colormap jet
        end 
        
        function draw_box(obj)
        
            v1 = [obj.grid_lower, obj.grid_lower];
            v2 = [obj.grid_lower, obj.grid_upper];
            v3 = [obj.grid_upper, obj.grid_upper];
            v4 = [obj.grid_upper, obj.grid_lower];
            box = [v1; v2; v3; v4; v1];
            plot(box(:,1), box(:,2), 'k', 'LineWidth', 4);
            
            
        end
        
       
        
        
        function arc(obj)
            comp.plotter.draw_arc(obj.max_r, obj.min_rad, obj.max_rad);
        end
        
     
        
        
        
        
        function print_info(obj)
            %PRINT_INFO Print the function associated with this plotter
        end
        
         
    end
      
end

