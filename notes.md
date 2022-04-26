

We might have 4 read-infiles

**species.yaml**

Stuff that almost never changes

- Species and atoms in the atmosphere, and their thermodynamic properties.

**settings.yaml**

Stuff that might change

- Planet parameters (gravity, etc.). What do to with water.
- opacities (k coeffs, rayleigh, ...)
- boundary conditions (species types, and BC types)

**star.txt**

- ascii file of solar flux

**atmosphere.txt**

- ascii file of initial atmospheric composition, and temperature.



# convection

Convection schemes are not time-dependent. They instantaneously transfer energy between layers.

We can add to the PDE describing temperature change, and put in quasi-realistic convection.

# notes on opacities

Clima uses lots of cross sections. For example, H2-H2 CIA is cross section I think. Also Rayleigh scattering is a cross section.

Also, Clima computes and interpolates this cross sections (and tau's) way more times than is necessary. It does it inside of the correlated-k loops. These calculations can happen before correlated-k loops.

To allow easier coupling with photochem, it might make sense to have the following.

One type of line opacity approximation:
- kdistributions

Different types of assumed continuum opacities:
- CIA
- Rayleigh
- photolysis xsections
- absorption xsections

For climate modeling only, there is no difference between photolysis and absorption cross sections. But distinguishing will make things easier for integration with photochemical modeling.

0-D, 1-D data will have two file types: ".txt" and ".dat". We will assume .txt is a text file, and that .dat is a fortran binary file. .txt will use photochem convention. We will make new convention for .dat.

2-D (xs as a function of lambda, T, P) data will only be binary: ".dat". Too confusing!


# Species types

Three different species types
- background filler gas
- fixed mixing ratio (fixed as a function of altitude)
- condensing gas (has boundary conditions)

default is fixed mixing ratio

Condensing gas + Temperature evolves through time.

Next the goal is to compute the heating rates. Given T, P, mixing ratio profile

Question: What are the independent and dependent variables? Input masses of H2O and CO2?

I think for now I will just consider only **fixed mixing ratios, with a background?**

P/(n*k) = T

We fix P at surface. We compute P above surface (according to T). We also compute T.

Only T evolves with time. This will work well for coupling to photochemistry.

The problem with this is that it isn't useful for all cases. Sometimes surface pressure is not an "input". Sometimes, we aim to compute surface pressure as H2O accumulates in the atmosphere. Post-impact atmospheres are a case in point.

Changing mixing ratio profiles in the climate model is challenging when coupling with photochemistry. Because chemistry will no longer influence the H2O profile in the upper atmosphere.

The solution is the track number densities in the photochemical model. This will permit the atmosphere to grow and shrink.

Let's start with fixed mixing ratios I think. It is easy, and still has big applications. An iterative approach can still be used to guesstimate the surface P. Basically, you can try to set surface P, which gives water at saturation, and the proper dry column abundances.




