clear;clc

wingroot=avlsect([0,0,0],0.2,0,'ag13');
wingtip=avlsect([0.05,0.5,0.1],0.1,0,'ag13');
wing=avlsurf("wing",2,[wingroot wingtip]);
tailroot=avlsect([1,0,0.2],0.1,0,'NACA0012');
tailtip=avlsect([1,0.25,0.2],0.1,0,'NACA0012');
tail=avlsurf("tail",-8,[tailroot tailtip]);

v=avlplane([wing tail]);

[CL CD Cl Cm Cn]=avl(v,operatingpoint(20,5,0))