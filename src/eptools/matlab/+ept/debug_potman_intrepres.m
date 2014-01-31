clear('pman');
pman.name = 'Laplace';
pman.size = 100;
pman.pars{1}=[0];
pman.pars{2}=[1.2];
% Comment out check at the end:
irp = ept.potman_intrepres(pman);
% Now run check
s = eptools_potmanager_isvalid(irp.potids,irp.numpot,irp.parvec, ...
			       irp.parshrd);

clear('pman2');
pman2{1} = pman;
pman2{2}.name = 'Gaussian';
pman2{2}.size = 200;
pman2{2}.pars{1}=randn(200,1);
pman2{2}.pars{2}=[1.5];
% Comment out check at the end:
irp2 = ept.potman_intrepres(pman2);
% Now run check
s2 = eptools_potmanager_isvalid(irp2.potids,irp2.numpot,irp2.parvec, ...
				irp2.parshrd);
