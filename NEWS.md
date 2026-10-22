# seqTrun 0.2.0

* Add simple and modified pseudo-observation regression for Cox models:
  `seqTrun.POreg.cox()` and `seqTrun.modPOreg.cox()`.
* Add simple and modified pseudo-observation regression for AFT models:
  `seqTrun.POreg.aft()` and `seqTrun.modPOreg.aft()`.
* Add the `geepack` dependency, help pages, and a simulated example.
* Replace placeholder definitions and global-data references in the supplied
  regression scripts with package functions. Use subject identifiers in GEE
  fits and keep regression helpers separate from the distribution estimators.
