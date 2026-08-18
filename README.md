# dark-star-evolution

The code models the coupled evolution of a dark star embedded within a baryonic star. The baryonic star is evolved with [MESA (https://docs.mesastar.org/), while the dark star is evolved using the instantaneous stellar structure of the baryonic star. The resulting dark star mass distribution determines the additional gravitational field acting on the baryonic star, which is then fed back into MESA.

```mermaid
flowchart TB
    MESA["MESA<br/>Baryonic Star"]

    EXTRACT["Extract Stellar Structure<br/>ρ(r), M(r), r"]

    DS["Dark-Star Evolution"]

    STRUCT["Solve Dark-Star<br/>Structural Equations<br/>dP/dr, dM/dr"]

    GRAV["Compute Gravitational Feedback<br/>g_extra(r)"]

    MESA --> EXTRACT
    EXTRACT -->|"Stellar environment"| DS
    DS --> STRUCT
    STRUCT -->|"Dark-star mass distribution<br/>M_DS(r)"| GRAV
    GRAV -->|"Additional gravitational acceleration<br/>on baryonic star"| MESA

    MESA -.->|"Next timestep"| EXTRACT
