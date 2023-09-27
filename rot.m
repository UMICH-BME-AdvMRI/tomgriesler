function R = rot(alpha, phase)
Rz = zrot(-phase);
Rx = xrot(alpha);
R = inv(Rz) * Rx * Rz;