classdef plotter
    %PLOTTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f              % f(z)        
    end
    
    properties(Hidden)
        C               % The color values used for contour projections and 2d color plots
        Z               % The matrix of input values for f(z)
        foZ             % The matrix of f(z)
        funName
        bounds
    end
    
    methods        
        function obj = plotter(anonFun, bounds, plotOnInstantiation)
            %PLOTTER Construct an instance of this class
            %   Pass an anonymous complex function to the constructor that
            %   you wish to be plotted
            assert(isa(anonFun, 'function_handle'), "First parameter must be a function handle")
            if nargin == 1
                bounds = [-10 10 -10 10];
                plotOnInstantiation = 1;
            elseif nargin == 2
                plotOnInstantiation = 1;
            end
            
            
            obj.f = anonFun;                        
            obj.bounds = bounds; % Set the bounds and compute f(z)

            obj = obj.updateName;

            if(plotOnInstantiation)
                obj.plot3;
            end
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
            figure('NumberTitle', 'off', 'Name', join(["3d plots of" obj.funName]));
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
        
        function mesh(obj) 
            %MESH_INTERSECT Plot a mesh of the real and imag parts of f(z)
            figure('NumberTitle', 'off', 'Name', "Mesh");
            r_mesh = mesh(real(obj.Z), imag(obj.Z), real(obj.foZ));
            hold on
            i_mesh = mesh(real(obj.Z), imag(obj.Z), imag(obj.foZ));
            hold off       
            axis square 
            my_map = [1 1 1; 0 0 0];
            colormap(gca, my_map);
            set(r_mesh, 'edgecolor', 'none', 'facecolor', [1 1 1]);
            set(i_mesh, 'edgecolor', 'none', 'facecolor', [0 0 0]);
            set(gca, 'Color', [.5 .5 .5]);
            xlabel('Re$(z)$', 'Interpreter', 'latex');
            ylabel('Im$(z)$', 'Interpreter', 'latex');
            zlabel('$f(z)$', 'Interpreter', 'latex');
            title("f(z)", 'Interpreter', 'latex');
        end

        function surf(obj) 
            %MESH_INTERSECT Plot a mesh of the real and imag parts of f(z)
            % figure('NumberTitle', 'off', 'Name', "Surf");
            r_surf = surf(real(obj.Z), imag(obj.Z), real(obj.foZ));
            hold on
            i_surf = surf(real(obj.Z), imag(obj.Z), imag(obj.foZ));
            hold off       
            axis square 
            my_map = [1 1 1; 0 0 0];
            colormap(gca, my_map);
            set(r_surf, 'edgecolor', 'none', 'facecolor', [1 1 1]);
            set(i_surf, 'edgecolor', 'none', 'facecolor', [0 0 0]);
            set(gca, 'Color', [.5 .5 .5]);
            xlabel('Re$(z)$', 'Interpreter', 'latex');
            ylabel('Im$(z)$', 'Interpreter', 'latex');
            zlabel('$f(z)$', 'Interpreter', 'latex');
            title(join(["$f(z) = " obj.funName "$"]), 'Interpreter', 'latex');
            legend('Re(x)', 'Im(y)', 'Interpreter', 'latex', 'Location', 'northeast')
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


        function obj = zoom(obj, xy_value) 

            obj.bounds = [-xy_value xy_value -xy_value xy_value];
            obj.plot3;

        end

        function L = limit(obj, z, r) 
            % I have to approach z from all directions...
            % I think that the best way 
            eps = 1e-5;
            if nargin == 2
                r = 1e-12;
            end
            theta = linspace(-pi, pi, 100);
            

            disk = r*exp(theta*i) + z;
            f_disk = arrayfun(obj.f, disk);
            disk_minus_d1 = f_disk - f_disk(1); % Subtract the first element, check if the resulting matrices entries are less then epsilon
            tiledlayout(1, 2);
            nexttile
            plot(disk);
            nexttile
            plot(f_disk);
            condition = disk_minus_d1 < eps;
            if all(condition) 
                L = f_disk(1);
            else 
                L = 'Limit not found';
            end
        end

        function z = zeros(obj)

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
            obj = obj.updateName;
            % disp("set.f called");
        end

        function obj = set.bounds(obj, bounds)
            obj.bounds = bounds;
            obj = obj.makeGrid;
            obj = obj.applyFun;
        end

        function obj = updateName(obj)
            fun_name = string(char(obj.f));
            fun_name = erase(fun_name, "@(z)");
            obj.funName = fun_name;
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

        function iterates(compFun, n)
        end


                 
    end

    methods(Static)

        function obj = trig
            b = 10;
            
            tiledlayout(2, 3);

            nexttile;
            obj = comp.plotter(@(z)sin(z), [-b b -b b], 0);
            obj.surf;

            nexttile;
            obj.f = @(z)cos(z);
            obj.surf;

            nexttile;
            obj.f = @(z)tan(z);
            obj.surf;

            nexttile;
            obj.f = @(z)asin(z);
            obj.surf;

            nexttile;
            obj.f = @(z)acos(z);
            obj.surf;

            nexttile;
            obj.f = @(z)atan(z);
            obj.surf;
        end

        function obj = poly
            b = 10;
            
            tiledlayout(2, 3);

            nexttile;
            obj = comp.plotter(@(z) z, [-b b -b b], 0);
            obj.surf;

            nexttile;
            obj.f = @(z)z^2;
            obj.surf;

            nexttile;
            obj.f = @(z)z^3;
            obj.surf;

            nexttile;
            obj.f = @(z)z^4;
            obj.surf;

            nexttile;
            obj.f = @(z)z^5;
            obj.surf;

            nexttile;
            obj.f = @(z)z^6;
            obj.surf;
        end

        function Im = julia(c, sb, res, max_iter, shouldPlot)
            %JULIA Plot the julia set for polynomial z^2 + c
            if nargin < 5
                shouldPlot = 1;
                if nargin < 4
                    max_iter = 250;
                    if nargin < 3
                        res = 1000;
                        if nargin < 2
                            sb = 2;
                            if nargin < 1
                                c = -0.8 + 0.156i;
                            end
                        end
                    end
                end
            end

            

            X = linspace(-sb, sb, res);
            Y = linspace(-sb, sb, res);
            [X, Y] = meshgrid(X, Y);
            Z = X + Y*i;
            C = ones(size(Z));

            r = (1 + sqrt(1 + 4*abs(c)))/2;

            f = @(z) (z^2 + c);

            for ii = 1:res
                for jj = 1:res
                    iter = 1;

                    z = Z(ii, jj);

                    for kk = 1:max_iter
                        if abs(z) > r
                            C(ii, jj) = kk/max_iter;
                            break
                        end
                        z = f(z);
                    end

                    if kk == max_iter
                        C(ii, jj) = 0;
                    end                    
                end
            end
            
            if (shouldPlot) 
                    
                disp("Julia set made, plotting")
                p = pcolor(real(Z), imag(Z), C);
                colormap(hot(256))
                set(p, 'edgecolor', 'none');
                colorbar;
                axis square
            end

            Im = C;
        end

        function juliaToGIF(nSteps, duration)
            if nargin < 2
                duration = 10;
                if nargin < 1
                    nSteps = 100;
                end
            end

            delay = duration/nSteps;
            theta = linspace(0, 2*pi, nSteps);
            C = exp(theta*i);
            map = hot(256);

            filename = 'test_julia.gif';

            for ii = 1:nSteps                
                Im  = comp.plotter.julia(C(ii), 1.8, 600, 150, 0);
                Im = uint8(rescale(Im, 0, 256));

                % filename = strcat("julia ", num2str(ii), ".png");
                % imwrite(Im, map, filename);
                if ii == 1
                    imwrite(Im, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.01);
                else
                    imwrite(Im, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
                end
            end

        end



    end





      
end

