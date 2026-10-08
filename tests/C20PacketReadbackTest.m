classdef C20PacketReadbackTest < matlab.unittest.TestCase
    properties
        Root
        Archive
        Report
    end
    methods (TestMethodSetup)
        function setupPacket(t)
            project=fileparts(fileparts(mfilename('fullpath')));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(project,'case_studies','bisb2026','scripts')));
            fixture=t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            t.Root=fixture.Folder; payload=fullfile(t.Root,'payload'); mkdir(payload);
            c20_v4_io('text',fullfile(payload,'sample.txt'),'synthetic readback only');
            writetable(c20_v4_io('inventory',payload),fullfile(payload,'FILE_MANIFEST.csv'));
            t.Archive=fullfile(t.Root,'packet.zip'); zip(t.Archive,{'*'},payload);
            t.Report=fullfile(t.Root,'report.json');
        end
    end
    methods (Test)
        function successPersistsReportAndCleans(t)
            result=c20_packet_readback(t.Archive,t.Report,@(p)struct('content',fileread(fullfile(p,'sample.txt'))));
            report=jsondecode(fileread(t.Report));
            t.verifyEqual(result.content,'synthetic readback only');
            t.verifyEqual(report.file_count,2);
            t.verifyEqual(report.zip_sha256,c20_v4_io('hash',t.Archive));
            t.verifyEqual(report.cleanup_status,'removed');
            t.verifyFalse(isfolder(report.temporary_directory));
        end
        function existingDirectoryIsUntouched(t)
            existing=fullfile(t.Root,'packet_readback'); mkdir(existing);
            c20_v4_io('text',fullfile(existing,'sentinel'),'keep');
            c20_packet_readback(t.Archive,t.Report,@(~)struct('passed',true));
            t.verifyEqual(fileread(fullfile(existing,'sentinel')),'keep');
        end
        function reportRefusesOverwrite(t)
            c20_v4_io('text',t.Report,'existing');
            t.verifyError(@()c20_packet_readback(t.Archive,t.Report,@(~)true),'c20:ReportExists');
            t.verifyEqual(fileread(t.Report),'existing');
        end
        function hashMismatchRetainsDiagnostics(t)
            payload=fullfile(t.Root,'payload'); c20_v4_io('text',fullfile(payload,'sample.txt'),'changed');
            delete(t.Archive); zip(t.Archive,{'*'},payload);
            t.verifyError(@()c20_packet_readback(t.Archive,t.Report,@(~)true),'c20:HashMismatch');
            t.verifyFalse(isfile(t.Report));
        end
        function callbackFailureRetainsDiagnostics(t)
            t.verifyError(@()c20_packet_readback(t.Archive,t.Report,@C20PacketReadbackTest.failVerification),'c20:TestFailure');
            t.verifyFalse(isfile(t.Report));
        end
        function traversalIsRejected(t)
            C20PacketReadbackTest.writeTraversalZip(t.Archive);
            t.verifyError(@()c20_packet_readback(t.Archive,t.Report,@(~)true),'c20:PathEscape');
            t.verifyFalse(isfile(fullfile(t.Root,'outside.txt')));
        end
        function manifestTraversalIsRejected(t)
            manifest=fullfile(t.Root,'payload','FILE_MANIFEST.csv');
            m=readtable(manifest,TextType='string'); m.path(1)="../outside.txt"; writetable(m,manifest);
            delete(t.Archive); zip(t.Archive,{'*'},fullfile(t.Root,'payload'));
            t.verifyError(@()c20_packet_readback(t.Archive,t.Report,@(~)true),'c20:FileSet');
            t.verifyFalse(isfile(fullfile(t.Root,'outside.txt')));
        end
    end
    methods (Static)
        function writeTraversalZip(path)
            stream=java.util.zip.ZipOutputStream(java.io.FileOutputStream(path));
            close=onCleanup(@()stream.close());
            stream.putNextEntry(java.util.zip.ZipEntry('../outside.txt'));
            stream.write(int8(42)); stream.closeEntry();
        end
        function result=failVerification(~)
            error('c20:TestFailure','Synthetic verifier failure');
            result=[]; %#ok<UNRCH>
        end
    end
end
