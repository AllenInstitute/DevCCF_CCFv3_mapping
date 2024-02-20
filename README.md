# DevCCF to CCFv3 mapping

Code for generating analyses for Figure 5 in Kronman et al. (doi: 10.1101/2023.09.14.557789)

### Requirements

#### System requirements

Mapping requires a single computer workstation with enough memory to load data.

#### Dependencies

All code is written in the python (3.11) environment

Required packages:

```
anndata >= 0.9 
anytree >= 2.7
plotly >= 5.10
seaborn >= 0.12 
SimpleITK >= 2.0
```

### Data

The `/data` subdirectory contains a set of .csv and .json files used to generate the figures. Users will need to acquire image volumes used in the manuscript, which can be found [here](https://kimlab.io/brain-map/DevCCF/)

### Level of support

We are not currently supporting this code, but simply releasing it to the community AS IS but are not able to provide any guarantees of support. The community is welcome to submit issues, but you should not expect an active response.
