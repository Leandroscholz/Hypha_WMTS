# Well-mixed tanks in series model for the growth of Aerial reproductive hyphae 

## Reference: W. Balmant, M. H. Sugai-Guérios, J. H. Coradin, N. Krieger, A. Furigo Junior, and D. A. Mitchell, “A Model for Growth of a Single Fungal Hypha Based on Well-Mixed Tanks in Series: Simulation of Nutrient and Vesicle Transport in Aerial Reproductive Hyphae,” PLoS One, vol. 10, no. 3, p. e0120307, Mar. 2015.

```main``` is the simplified model implementation in which the simulation begins with N reactors (N > Nv) and the tip-tank reactor extends indefinitely until the defined time span ```tspan```

```mainWithEvent.m``` is the simplified model implementation in which the simulation begins with N reactors (N > Nv) and the simulation ends on time ```t``` in which the tip-tank reactor reaches twice its initial length. 

```mainWithEventAndLoop.m``` is the complete model implementation. N can be any positive integer and the simulation runs until the maximum number of reactors ```NMax``` (defined in row 19 as input parameter) is reached.
