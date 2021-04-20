function mribx = SetBox( segmentedmri, bx )

Tr = eye(4);
Tr(1,4) = bx(1,1)-1;  Tr(2,4) = bx(2,1)-1;  Tr(3,4) = bx(3,1)-1;  

mri_sgbx = segmentedmri.anatomy(bx(1,1):bx(1,2), bx(2,1):bx(2,2), bx(3,1):bx(3,2) );

mribx = [];
mribx.dim       = size(mri_sgbx);
mribx.transform = Tr;
mribx.segTr     = segmentedmri.transform;
% mribx.coordsys  = segmentedmri.coordsys;
mribx.unit      = segmentedmri.unit;
mribx.anatomy   = mri_sgbx;
mribx.box       = bx;

% ft_sourceplot([], mribx)


