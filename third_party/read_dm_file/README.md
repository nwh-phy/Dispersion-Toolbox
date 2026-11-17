# ReadDMFile

Fred Sigworth, version 1.1.0 (2013), BSD license in LICENSE.txt.
Source: https://www.mathworks.com/matlabcentral/fileexchange/43005-read-dm3-and-dm4-image-files
Downloaded 2026-11-17 from version 2 ZIP.

Local changes: optional third argument selects the ImageList member (default 2
preserved); support fourth dimension; fourth output exposes selected calibration
values (8: energy origin; 9-11: second-axis scale, units, origin).
The Bi slot scan is member 1 and dimensions E x detector channel x scan-x x scan-y.
This upstream reader assumes little-endian data and machine; tested here only on
ARM macOS with the little-endian Nion-exported DM4. Other DM encodings unverified.
