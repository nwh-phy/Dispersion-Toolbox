classdef qe_plot_helpers
%QE_PLOT_HELPERS  Reusable plotting utilities for q-E browser.
%
%   Static methods for branch scatter/errorbar plotting, dispersion
%   fit overlay, map display preparation, and shared color palettes.
%
%   Usage:
%     qe_plot_helpers.plot_branch_scatter(ax, br, col, label)
%     colors = qe_plot_helpers.branch_colors()
%
%   See also: interactive_qe_browser, qe_gamma_dashboard

    methods (Static)

        function colors = branch_colors()
            %BRANCH_COLORS  Standard 5-color palette for branches.
            colors = [0.1 0.5 0.9; ...
                      0.9 0.3 0.1; ...
                      0.2 0.7 0.3; ...
                      0.6 0.2 0.8; ...
                      0.9 0.6 0.1];
        end


        function col = branch_color(b)
            %BRANCH_COLOR  Get color for branch index b (wraps around).
            colors = qe_plot_helpers.branch_colors();
            col = colors(mod(b-1, size(colors,1))+1, :);
        end


        function plot_branch_scatter(ax, br, col, label)
            %PLOT_BRANCH_SCATTER  Plot branch points with optional CI error bars.
            %
            %   Plots q vs E (columns 1 and 2 of br). If columns 6-7 contain
            %   valid CI bounds, error bars are drawn. Otherwise, filled scatter.
            if size(br, 2) >= 7 && ~all(isnan(br(:,6)))
                E_err = (br(:,7) - br(:,6)) / 2;
                E_err(isnan(E_err)) = 0;
                errorbar(ax, br(:,1), br(:,2), E_err, ...
                    'o', 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 4, 'LineWidth', 0.8, 'CapSize', 3, ...
                    'DisplayName', label);
            else
                scatter(ax, br(:,1), br(:,2), 30, col, 'filled', ...
                    'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', label);
            end
        end


        function plot_fit_curve(ax, fit_result, col, branch_idx)
            %PLOT_FIT_CURVE  Overlay a dispersion fit curve on axes.
            if isempty(fit_result), return; end

            if isfield(fit_result, 'model_label')
                lbl = sprintf('Fit B%d: %s  R²=%.3f', ...
                    branch_idx, fit_result.model_label, fit_result.R_squared);
            elseif isfield(fit_result, 'rho0')
                lbl = sprintf('Fit B%d: ρ₀=%.1fÅ  E_f=%.0fmeV', ...
                    branch_idx, fit_result.rho0, fit_result.E_flat_meV);
            else
                lbl = sprintf('Fit B%d: R²=%.3f', branch_idx, fit_result.R_squared);
            end

            plot(ax, fit_result.q_fit, fit_result.E_fit, '-', ...
                'Color', col, 'LineWidth', 2.5, ...
                'DisplayName', lbl);
        end


        function [display_map, display_clim] = prepare_map_display(map_values, mode_name, display_style)
            %PREPARE_MAP_DISPLAY  Scale map values for imagesc display.
            %
            %   [display_map, clim] = prepare_map_display(values, mode, style)
            %
            %   mode_name    — 'display' or other trace mode
            %   display_style — 'physical' or 'normalized'
            if nargin < 2, mode_name = "display"; end
            if nargin < 3, display_style = "physical"; end

            display_map = double(map_values);
            finite_values = display_map(isfinite(display_map));
            if isempty(finite_values)
                display_clim = [0 1];
                display_map(:) = 0;
                return
            end

            switch mode_name
                case "display"
                    if strcmpi(display_style, "normalized")
                        lower_limit = max(min(finite_values), 0);
                        upper_limit = prctile(finite_values, 99.8);
                        if ~isfinite(upper_limit) || upper_limit <= lower_limit
                            upper_limit = max(finite_values);
                        end
                        if ~isfinite(upper_limit) || upper_limit <= lower_limit
                            upper_limit = lower_limit + 1;
                        end
                        display_map = max(display_map, lower_limit);
                        display_clim = [lower_limit, upper_limit];
                    else
                        lower_limit = prctile(finite_values, 0.5);
                        upper_limit = prctile(finite_values, 99.8);
                        if ~isfinite(lower_limit) || ~isfinite(upper_limit) || upper_limit <= lower_limit
                            lower_limit = min(finite_values);
                            upper_limit = max(finite_values);
                        end
                        scale = max(upper_limit - lower_limit, eps);
                        clipped = max(display_map - lower_limit, 0);
                        stretch = 8;
                        display_map = asinh(stretch * clipped / scale) / asinh(stretch);
                        display_clim = [0 1];
                    end

                otherwise
                    center_value = median(finite_values, "omitnan");
                    spread = prctile(abs(finite_values - center_value), 99);
                    if ~isfinite(spread) || spread <= eps
                        spread = max(abs(finite_values - center_value));
                    end
                    if ~isfinite(spread) || spread <= eps
                        spread = 1;
                    end
                    stretch = 4;
                    display_map = asinh(stretch * (display_map - center_value) / spread) / asinh(stretch);
                    display_clim = [-1 1];
            end
        end

    end
end
