# bHIVE

## Artifical Immune Network-based Machine Learning Pacakge

<img align="right" src="https://github.com/ncborcherding/bHive/blob/main/www/bhive_hex.png" width="305" height="352">

### Introduction

Artificial Immune Systems (AIS) in general mimic mechanisms found in the natural immune system. Examples include:

* Clonal Selection: B-cells (or T-cells) that bind an antigen with high affinity proliferate and mutate (affinity maturation).
* Immune Network Theory: Jerne’s idiotypic network hypothesis states that antibodies can recognize not only antigens but also other antibodies, giving rise to a “network” of interactions.

An artificial immune network (AIN) is typically an unsupervised or semi-supervised algorithm that simulates how these immunological elements (antibodies and their interactions) evolve, cluster data, and detect “novel” patterns.

**Key Immune-Inspired Operations**

* Affinity – A measure of how well an antibody (candidate solution) “recognizes” or “matches” a pattern (e.g., input data). Often a similarity or distance metric.
* Clonal Selection and Expansion – Antibodies that exhibit high affinity to patterns get cloned and undergo mutation (exploration around the “parent” solution).
* Network Suppression – Redundant or overly similar antibodies are suppressed (pruned), keeping the set of antibodies diverse.
