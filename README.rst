.. role::  raw-html(raw)
    :format: html

*JEVOL* : Angular moment evolution code
=======================================

This numerical code is developed to follow the evolution of the rotation rate of low mass stars 0.5 < M\ :sub:`star`\/M\ :sub:`sun`\  < 1.0 from the PMS to the end of the MS.

*JEVOL* is writen in Fortran77 

The input parameters are stored in the rotevoldec.par file

Input parameters: 

1) pinit: the initial rotation period (P\ :sub:`rot`\) in days
2) mass: the mass of the star
3) stellar model : the stellar evolution model used
4) ksk, kmm, ksc, kmp: Skumanich, Mayor-Mermilliod, Schatzmann, and  Matt braking law constant constants
5) K, K1MP, K2MP, m: Matt et al. (2012) braking law constants
6) taudec: core-envelope coupling timescale :raw-html:`&tau;`:sub:`c-e` in year
7) tdisk: star-disk interaction timescale :raw-html:`&tau;`:sub:`disk` in year 
8) brklaw: selection the braking law (0 = Matt et al. 2012); 1-2 = Reville et al. (2015); 3 = MAtt et al. (2015)) 
9) K3, K4, mvic: Reville et al. (2015) braking law constants 


.. code-block:: bash

    make rotevoldec
    ./rotevoldec.e

