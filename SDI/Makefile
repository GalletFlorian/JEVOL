##############rotevoldec
rotevoldec.o : rotevoldec.f
	gfortran rotevoldec.f -O -c 

magnbrak.o : magnbrak.f
	gfortran magnbrak.f -O -c 

magnbrakcl.o : magnbrakcl.f
	gfortran magnbrakcl.f -O -c 

mdot.o : mdot.f
	gfortran mdot.f -O -c 

interspline.o : interspline.f
	gfortran interspline.f -O -c 

spline.o : spline.f
	gfortran spline.f -O -c 

splint.o : splint.f
	gfortran splint.f -O -c 

interlin.o : interlin.f
	gfortran interlin.f -O -c 

#ROTEVOLDEC = rotevoldec.o magnbrak.o interspline.o spline.o splint.o
ROTEVOLDEC = rotevoldec.o magnbrak.o magnbrakcl.o interlin.o mdot.o 
rotevoldec : $(ROTEVOLDEC)
	 gfortran -o rotevoldec.e $(ROTEVOLDEC)

