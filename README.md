# TCBFlash
This project aims to measure the environmental quenching timescale for
the classical Local Group dwarf galaxies through comparison to the
ELVIS suite of dark matter only simulations
## Radial Model
This quenching model uses a fixed radial distance as a function of the
virial radius, inside which, all subhalos are labeled as
quenched. This quenching radius can be varied so as to match the
observed quenched fraction in the Local Group.
## Time Model
This quenching model uses a fixed quenching timescale for all of the
subhalos in the ELVIS simulation. Such that any subhalo will be
labeled as quenched once it has spent that amount of time inside the
virial radius.
* The clock starts at virial infall
* A subhalo "quenched" when the clock reaches the quenching timescale.
* The quenching timescale is varied so as to match the observed
  quenched fraction in the Local Group.
