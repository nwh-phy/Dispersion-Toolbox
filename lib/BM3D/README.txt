BM3D — Block-Matching and 3D Filtering
=======================================

Original algorithm by K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian.
Tampere University of Technology, Finland.
http://www.cs.tut.fi/~foi/GCF-BM3D

BM3D_QRS wrapper modified by Ruishi Qi, Peking University.

This code is redistributed for academic/research use only.
If you use BM3D in your work, please cite:

  K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian,
  "Image denoising by sparse 3D transform-domain collaborative filtering,"
  IEEE Trans. Image Process., vol. 16, no. 8, pp. 2080-2095, Aug. 2007.

Files included:
  BM3D.m          — Core BM3D algorithm
  BM3D_QRS.m      — Automatic noise estimation wrapper
  bm3d_thr.*      — Block-matching hard thresholding (compiled MEX)
  bm3d_wiener.*   — Wiener collaborative filtering (compiled MEX)

Platform MEX binaries:
  .mexw64    — Windows 64-bit
  .mexa64    — Linux 64-bit
  .mexmaci64 — macOS 64-bit (Intel)
