classdef Run590GuiHistoryAreaOptionsTest < matlab.unittest.TestCase
    methods (Test)
        function acceptsQRangeAndOutputTagOptions(testCase)
            caseRoot = fileparts(fileparts(mfilename('fullpath')));
            scriptPath = fullfile(caseRoot, 'scripts', ...
                'run_590_gui_history_area_analysis.m');

            testCase.verifyTrue(isfile(scriptPath));
            src = fileread(scriptPath);

            testCase.verifyNotEmpty(strfind(src, ...
                'options.qRangeOverride_Ainv'));
            testCase.verifyNotEmpty(strfind(src, ...
                'options.outputTagSuffix'));
            testCase.verifyNotEmpty(strfind(src, ...
                'options.peakModelOverride'));
            testCase.verifyNotEmpty(strfind(src, ...
                'local_apply_run_options_to_session'));
            testCase.verifyNotEmpty(strfind(src, ...
                'snap.qStart = options.qRangeOverride_Ainv(1);'));
            testCase.verifyNotEmpty(strfind(src, ...
                'snap.qEnd = options.qRangeOverride_Ainv(2);'));
            testCase.verifyNotEmpty(strfind(src, ...
                'snap.peakModel = char(string(options.peakModelOverride));'));
        end

        function defaultsQRangeToQ015WhenUnspecified(testCase)
            caseRoot = fileparts(fileparts(mfilename('fullpath')));
            scriptPath = fullfile(caseRoot, 'scripts', ...
                'run_590_gui_history_area_analysis.m');

            testCase.verifyTrue(isfile(scriptPath));
            src = fileread(scriptPath);

            testCase.verifyNotEmpty(strfind(src, ...
                'options.qRangeOverride_Ainv (1,2) double = [-0.15 0.15]'));
            testCase.verifyNotEmpty(strfind(src, ...
                'options.peakModelOverride {mustBeTextScalar} = ""'));
            testCase.verifyNotEmpty(strfind(src, ...
                'if all(isfinite(options.qRangeOverride_Ainv))'));
            testCase.verifyNotEmpty(strfind(src, ...
                'session.output_tag = [session.output_tag char(options.outputTagSuffix)];'));
        end
    end
end
