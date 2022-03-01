

We might have 3 read-infiles

**species.yaml**

- Species and atoms in the atmosphere, and their thermodynamic properties.

**settings.yaml**

- Planet parameters (gravity, etc.)
- opacities (k coeffs, rayleigh, ...)

**star.txt**

- ascii file of solar flux


# convection

Convection schemes are not time-dependent. They instantaneously transfer energy between layers.

We can add to the PDE describing temperature change, and put in quasi-realistic convection.





# notes on opacities

Clima uses lots of cross sections. For example, H2-H2 CIA is cross section I think. Also Rayleigh scattering is a cross section.

Also, Clima computes and interpolates this cross sections (and tau's) way more times than is necessary. It does it inside of the correlated-k loops. These calculations can happen before correlaed-k loops.






