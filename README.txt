a la mode: mode matching and beam propagation solutions for MATLAB
nicolas smith

It requires MATLAB 2008a, you should add the alm folder to your path
(the @ directories should be in a folder that's in your path, but they
shouldn't be in the path), or do your work in that folder.

type 'help alm' in MATLAB to begin

example scripts (located in examples/):
pslmode.m - showing some basic uses
choosecompsexample.m - shows how to use the chooseComponents method
beamfitexample.m - uses the beam width fitting function
customcostexample.m - shows how to define a custom cost function when optimizing your beam path
