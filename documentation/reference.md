Yell Reference
=====
Refers to the program version Yell 1.0
Latest manual update: February 6th, 2014

[TOC]

Installation
=========

Yell is distributed as a self-contained executable. Just copy it in a folder where the shell can find it. Instructions for building the binary from the source files will be added later.
    
How to run Yell
===========

Yell is a simple console application. In order to run it, open a terminal, navigate to the folder with the `model.txt` file and type  

    yell.exe

(on Windows) or 

    yell 
    
(on Mac).

Input files
======

## Yell data format for arrays {#array-format}

All input arrays should be presented in  [hdf5](http://www.hdfgroup.org/HDF5/whatishdf5.html) formatted files. The hdf5 format is commonly used in scientific applications to hold large data arrays. It is not editable by hand, but the libraries for preparing and manipulating the files are available for most of the software platforms including Matlab, Python, Java, C++ and R.

Each Yell-formatted hdf5 files should have the following structure: the root `/` of the file is expected to have attribute 'filetype' equal to 'Yell 1.0' and the dataset should be located at the path '/data' within the file. The array should have [row-major order](http://en.wikipedia.org/wiki/Row-major_order). Though, in principle, hdf5 files may contain several datasets, Yell only expects one dataset per file.

Here is code to save file in this format in Matlab:

```
function write4yell(filename,dataset)
  dataset=permute(dataset,[3 2 1]); % Matlab uses column-major order, thus we need to permute dimensions
  h5create(filename,'/data',size(dataset));
  h5write(filename,'/data',double(dataset));
  h5writeatt(filename,'/','filetype','Yell 1.0');
```

## Allowed array input files 

The following files are used 

* **experiment.h5** File containing the experimental diffuse scattering (in reciprocal space). Required for refinement.
* **weights.h5** File containing the least-squares weights (in reciprocal space). Optional. If not present unit weights are used.
* **reciprocal_space_multiplier.h5** File containing the array that will be multiplied to the calculated diffuse scattering (in reciprocal space). Useful for zeroing-out the model intensities in places where experiment was not measured. Optional.
* **pdf_multiplier.h5** File containing the array that will be multiplied to the calculated PDF. Useful to simulate resolution function effects. Optional.

All arrays must have the same shape along all dimensions

The following indicates at which stages the arrays are used:  
1. Diffuse scattering is calculated.  
2. If available, PDF multipliers are applied in PDF space 
3. If available, reciprocal space multipliers are applied in reciprocal space 
4. The resulting dataset is called 'model intensities'  

During a refinement $weights*(I_{experiment}-I_{model})$ is used for minimization. 

No input arrays files are required if the refinement option is switched off (see below).


## model.txt

This file contains the definition of the local structure model. For description see [reference](#reference).

Output files
======


The output files are also in hdf5 format. The following shows a Matlab function which reads in the data:

```
function res=read_from_yell(filename)
    res = h5read(filename,'/data');
    res = permute(res,[3 2 1]); #Transform from row-major order
```


## Reciprocal space datasets
* **model.h5** Diffuse scattering calculated from the model. If a refinement is performed, `model.h5` is calculated from the final parameters; refined intensities are scaled to the experimental intensities.
* **average.h5** Full scattering calculated from average pairs. Only takes pairs touched by correlations as defined in model.txt, i.e. broad Bragg peaks and truncation ripples are possible. Scaled to the experimental intensities.
* **full.h5** Full scattering from real structure pairs. Only takes those pairs which were touched by correlations. Scaled to the experimental intensities.

## ∆PDF space datasets
All of this datasets are calculated only if [diffuse scattering grid](#DiffuseScatteringGrid) allows FFT. 
The datasets are scaled to be in $e^2/A^3$, if the [multiplicity](#Multiplicity) is set correctly.

* **delta-pdf.h5** Modeled ∆PDF.   
* **exp-delta-pdf.h5** Experimental ∆PDF, i.e. the normalized Fourier transform of the diffuse intensities provided by experiment.h5. 
* **delta-delta-pdf.h5** ∆$^2$PDF = ∆PDF$_{exp}$ - ∆PDF$_{model}$.  


Keywords reference {#reference}
========
## File format

The file `model.txt` is a plane text file which describes the crystal model. The file contains the following sections:

* **preamble** defines parameters which influence the calculation
* **UnitCell** section defines the average crystal structure
* **Modes** section defines the possible motions of the atoms or molecules
* **Correlations** section defines the short-range ordering in the crystal
* in **epiloge** the user may request Yell to print various expressions

The input file contains square brackets `[` `]` for grouping various elements, assignments for providing aliases to the various parts of the structure and also variables and arithmetic expressions in order to express constraints.

The parser in Yell is quite flexible considering white spaces. Space, tabulation, newline character and comment are all considered whitespace. In most parts of the input files the whitespace characters are not required, but user can add them to improve the readability.

## Comments

Comments in Yell start with a sharp sign `#` and go until the end of the line.

    #this is a comment
    Cell 1 1 1 90 90 90  #this is also comment

## Assignments

The parts of the sturucture - atoms, groups of atoms, variants and modes can be given aliases for the later reference in the model. The aliases are expressed as assignments

    Cu1_alias = Cu 1  0 0 0  0.01

## Formulas and variables

### Arithmetic expressions {#expressions}

The parser of Yell allows the arithmetic expressions in the place of almost any number in the input file. The expressions will be calculated during diffuse scattering calculation according to usual mathematical rules. Thus the following line

     Cell   6-1   1+2*2   4       9*(15-5)   3*3*2*5   120

is equivalent to 

    Cell  5 5 4    90 90 120

Expressions can contain the operators `+` `-` `*` `/`, brackets `()`, and special mathematical functions `exp` `log` `sin` `cos` `sqrt` `abs` `mod` and `pow`. The whitespace characters are **not** allowed within formulas and will invoke an error.

### Variables

User can define variables to use in arithmetic expressions. The definition of variables is similar to many programming languages:

    a=1;
    b=12*3;

Variable definition should end with semicolon `;`. The whitespace characters are **not** allowed in the variable definition, e.g. the expression `b = 12 * 3;` is illegal.

Valid names for variables start with upper or lower case letters, names may also contain numbers and underscore. The names are case sensitive. Here are the examples:

    A=1;
    shift_along_111=1;

The variables can then be used in any arithmetic expression. In any case *float* algorithms are used.

    a=12.1;
    c=a/2;
    alpha=90;
    gamma=60;

    Cell a a c alpha alpha gamma*2

Variables can be defined only in the "top level" of the input file, i.e. in  preamble, epilogue and at the top levels of the definitions of `UnitCell` `Modes` and `Correlations`. Here is an example which explicitly describes where the variables can and cannot be defined

    a=1;
    Cell 10 10 10 90 90 90
    b=1; #Here is ok
    DiffuseScatteringGrid -10 -10 -10  0.5 0.5 0.5  20 20 20
    UnitCell
    [
      c=1; #Here is ok
      Var = Variant
      [
        #but not here
        (p=0.5) #this is ok because it is a keyword expression and not an assignment to a variable 'p'
        [
          #not here
          Ag 1  0 0 0   0.02
        ]
        (p=0.5)
        Void
      ]
      d=1; #Here is ok
     ]
    Modes
    [
    e=1; #Here is ok
    ]
    Correlations
    [
      f=1; #Here is ok
      [(0,0,0)
        #not here
        SubstitutionalCorrelation(Var,Var,0.5)
      ]
    ]
    g=1; #Here is ok
    Print "g=" g

### Refinable variables
      
User can also define a variable as a *refinable variable*. Such variables will be used in the least squares refinement procedure. Their value will be optimized in order to fit experimental data. Other variables, which were not defined as refinable, but depend on refinable variables in the course of refinement will always be updated, based on the actual value of refinable variables they depend upon.

Refinable variables should be refined in the preamble of the file, with `RefinableVariables` keyword. Example:

    RefinableVariables
    [
    a=10; #This variable will be changed in the course of a refinement
    ]
    b=a*2; #This variable will change its value if the value of variable 'a' is modified during the refinment. 

Refinable variables can **not** be initialized by expressions:

    RefinableVariables
    [
    a=10/2; #Error
    ]

### Special functions that may be used in expressions
`exp(x)` $=e^x$  
`log(x)` Natural logarithm. Invokes error, if $x \leq 0$.  
`sin(x)` sine  
`cos(x)` cosine  
`sqrt(x)` $=\sqrt{x}$. Invokes error if $x<0$.  
`abs(x)` $=|x|$  
`mod(x,y)` $= x \mod y$  
`pow(x,y)`$=x^y$  

## Preamble

The first part of `model.txt` contains general settings of the calculation, like definition of the unit cell, grid for diffuse scattering calculation, point group, definition of refinable variables, the calculation method and several parameters which can affect the calculation and refinement process.

The mandatory fields are `Cell` `DiffuseScatteringGrid` and `PointGroup`, other keywords can be omitted. In such cases the program will use default values.


### Cell
Mandatory field.  
Defines the unit cell of the crystal.
Format: Cell $a$ $b$ $c$ $\alpha$ $\beta$ $\gamma$.  
Units: Ångströms and degrees.
Example:

    Cell 5.406 5.406 5.406  90 90 90

### DiffuseScatteringGrid {#DiffuseScatteringGrid}
Mandatory field.  
Defines the grid on which the diffuse scattering should be calculated. In refinements experimental data and weights are to be provided with the same grid definition. 

The format is the following:

    DiffuseScatteringGrid lower_limit_x lower_limit_y lower_limit_z  
                          step_size_x step_size_y step_size_z  
                          number_of_pixels_x number_of_pixels_y number_of_pixels_z

where `lower_limits` are the minimal `h` `k` and `l` indices of the grid, `step_sizes` are the distances between lattice points in units of `h` `k` and `l`, and `number_of_pixels` are the total number of pixels along each direction.

The main axes of the grid can only be defined along the main axes of the unit cell. If one needs to calculate the diffuse scattering grid along special directions, like, for example `111`, one has to transform the unit cell accordingly.

Yell uses a Fast Fourier Transform (FFT) algorithm for switching between reciprocal and PDF space. This adds two constraints on the diffuse scattering grid:

1. number of pixels along each dimension must be even
2. the origin of reciprocal space must be in the pixel with coordinates `(number_of_pixels_x/2+1, number_of_pixels_y/2+1, number_of_pixels_z/2+1)`

The conditions are a bit unusual for anyone who did not have experience with the FFT before. The arrays with uneven number of pixels and the central pixel in the center seems more natural. In order to use FFT such "natural" datasets should be stripped of the last planes along x y and z direcitons.

The FFT algorithm allows to calculate cross-sections through the center of reciprocal space. Such cross-sections in reciprocal space correspond to the projection of the whole structure to a plane, or a line in PDF space. If the cross-section is calculated the number of pixels is allowed to be equal to 1 along some axis, then the `step_size` along the corresponding dimension is ignored, the `lower_limit` should be equal to 0. 

Examples:

    DiffuseScatteringGrid -5 -5 -5  0.1 0.1 0.1  50 50 50   #Usual three-dimensional case
    DiffuseScatteringGrid -10 -10 0   0.1 0.1  1   20 20 1  #hk0 section of diffuse scattering

If the grid is not consistent with the above mentioned rules, diffuse scattering can be calculated using the [*direct*](#CalculationMethod) calculation algorithm and refinement can still be performed. However, the FFT algorithm is required if you want to

* apply **pdf_multipliers.h5**
* use [fast diffuse scattering calculation algorithm](#CalculationMethod)
* obtain **delta-pdf.h5**, **exp-delta-pdf.h5** and **delta-delta-pdf.h5**.

The following example only works with 'exact' calculation method:

    DiffuseScatteringGrid 0 0 0  0.2 0.2 0.2  30 30 30      #Only works with 'exact' calculation method

### LaueSymmetry
Mandatory field.  
Possible values: `m-3m` `m-3` `6/mmm` `6/m` `4/mmm` `4/m` `-3:R` `-3:H` `-3m:H` `-3m:R` `mmm` `2/m` `2/m:b` `-1`

Defines the Laue group of the crystal. 

Laue groups which have both rhombohedral and hexagonal settings are noted by the `:R` or `:H` letters. The group `2/m` also have two standard settings, with unique `c` and with unique `b`. They are called `2/m` and `2/m:b` respectedly. Non-crystallographic symmetry, e.g. five-fold rotation axes, are not supported by the program, but can be provided explicitly by the definition of correlations (see below).  

### Scale {#scale}
default value: 1

Scaling coefficient between model and experiment; always refined.

### RefinableParameters
default: empty `[ ]`

Registers variables in the square brackets to be used during refinement.  The variables cannot be initialized with expressions.

Example:

    RefinableVariables
    [
      a=1;
      b=12;
    ]

The scale coefficient, which is also refined, should not be defined here, since it has a special [keyword](#scale).

### RecalculateAverage
possible values: `true` `false`  
default value: `true`

Controls whether the average PDF of the structure should be recalculated during refinement. If set to `false` the average PDF is only calculated  in the beginning of refinement and kept unchanged throughout the refinement.

In cases, when non of the refined variables may influence the average interatomic vectors, because the average structure is very well known, the recalculation can be turned off. This speeds up the refinement by the factor of two. It is safe to do so, when refinable variables do not appear in sections `UnitCell`, `Modes`, and in `Correlations` in the brackets `(1,0,0)` and in `Multiplicity`.

### DumpPairs
possible values: `true` `false`  
default: `false`

When turned on, prints all the interatomic pairs used to calculate diffuse scattering in the program output. 

### PrintCovarianceMatrix
possible values: `true`, `false`
default: `false`

When turned on, prints the full refinement covariance matrix in the program output.

### CalculationMethod {#CalculationMethod}
possible values: `exact` `approximate`   
default value: `exact`

Selects the diffuse scattering calculation method. Yell currently has two calculation methods:

* `exact ` uses direct sum over all pairs, as shown in [this formula](#TheFormula). Calculation is done in reciprocal space.
* `approximate` goes through real space. This works much quicker because for each interatomic pair it only calculates a certain block where the signal is significant (controlled by flag [`FFTGridSize`](#FFTGridSize)) and ignores most of the PDF space , where the signal is almost zero. Details are described in ([Simonov et al. in. prep.](#YellPaper))

The approximate algorithm provides a significant speedup, but introduces errors. It is adviced to properly set up the `approximate` algorithm in the first stages of refinement, but always check the results with `exact` refinement before publishing. 

### Parameters controlling approximate diffuse scattering calculation

#### FFTGridSize {#FFTGridSize}
expects: three numbers  
default: `16 16 16`  

Controls the size of parallelepiped in pixel units to calculate the PDF signal of a pair. Bigger values give better accuracy, smaller provide more speed. The values should not be bigger than the size of the dataset along the same dimension.

#### FFTGridPadding
expects: three numbers  
default: `0 0 0`

Defines the padding in the approximate algorithm.

Prior to PDF calculation, Yell extends reciprocal space by the selected number of pixels; after the calculation the padding is cut.

Frequently, the FFT approximation method has the strongest errors close to the edges of reciprocal space. The padding expands the reciprocal space putting the errors outside the region of interest.

Details are described in ([Simonov et al. in. prep.](#YellPaper)).

#### PeriodicBoundaries {#PeriodicBoundaries}
expects: three booleans  
default: `true true true`

During the `approximate` diffuse scattering calculation allows not to periodically repeat the ∆PDF signals which lie outside the calculated PDF grid along specified directions.

**Why this might be interesting**  
If the correlations are very long in PDF space, one is forced to reconstruct diffuse scattering on a very fine grid, and use a huge array for ∆PDF refinement. However, sometimes the values of the correlations with long correlation vectors are not interesting. In such case there is a trick one can use to reduce the amount of memory Yell requires for refinement.

Experimental data is prepared in a special way. From the reconstruction on the fine grid one calculates the ∆PDF. Then cuts a central part of the ∆PDF, containing the region of interest. And finally back-transforms the part of the ∆PDF back to reciprocal space.

This procedure introduces strong truncation ripples in reciprocal space, but preserves the ∆PDF. The ∆PDF can then be refined in Yell with the FFT method. The obtained correlations will be accurate, though the diffuse scattering will be modeled poorly. 

In such case interatomic pairs which lie close to the edge of calculated ∆PDF region. The part of the ∆PDF signal from such pairs will lie inside the calculated region, and part of the signal will lie outside the calculated region. The *direct* calculation method and the *approximate* algorithm with default settings will wrap the ∆PDF signals which lie outside the calculated region into the calculated area. Turning periodic boundaries off along these dimensions avoids this problem.

Thus, whenever the above mentioned method of preparing diffuse scattering is used, the periodic boundaries should always be turned off along the dimensions along which the ∆PDF map was cut.

### ReportPairsOutsideCalculatedPDF
possible values: `true`, `false`
default: `false`

Reports all the pairs which lie outside calculated ∆PDF map. 

If an interatomic pair fall outside the calculated ∆PDF map, its signal will be periodically wrapped inside the ∆PDF map. For diffuse scattering comprising layers and streaks this is a valid behavior. However, if the diffuse scattering is not broad along some dimension, but certain pairs fall outside the calculated ∆PDF map along that dimension, this is an indication that diffuse scattering should be reconstructed on finer grid.

During the refinement with a [specially prepared dataset](#PeriodicBoundaries), with some [`PeriodicBoundaries`](#PeriodicBoundaries) turned off, the signal from most of the interatomic pairs which lie outside the calculated region will be silently discarded. Turning on this option allows to find unexpected behavior caused by this feature.

### Parameters controlling the least square refinement
#### Refine
Possible values: `true`, `false`  
Default: `true`

Specifies whether the program should refine diffuse scattering. If set to `false`, the program will just calculate diffuse scattering from a given model.

#### MaxNumberOfIterations
possible values: integer
default: 1000  

Defines the maximum number of iterations the refinement algorithm is allowed to run.

#### MinimizerTau
default: 1E-03

Defines the initial Levenberg-Marquandt damping parameter. Note that this value works differently, than the `damping` in standard crystallographic packages. For more information see this [wikipedia article](http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) and also somewhat minimalistic [manual](http://users.ics.forth.gr/~lourakis/levmar/) of the levmar library. 

#### MinimizerThresholds
possible values: three numbers
default: 1E-17 1E-17 1E-17

Defines the thresholds to detect convergence of least squares procedure. The thresholds are for the size of the gradient, minimization step size, the third criterion works if the experimental and diffuse scattering are the same (say, for synthetic data). For more info see `\epsilon1` `\epsilon2` and `\epsilon3` in [levmar documentation](http://users.ics.forth.gr/~lourakis/levmar/).

#### MinimizerDiff
possible values: one number  
default: 1E-06  

Yell calculates derivatives by finite difference method. This value provides the increment that is used for the calculation of derivatives.

## Unit cell

This section of the input file defines the contents of the crystal's average unit cell. 

In Yell, the average structure is defined in a hierarchical manner: the unit cell consists of topological sites, which are occupied by atoms or groups of atoms. The topological sites are called `Variants`and can be disordered (occupied with a certain probability by more than one group of atoms). Groups of atoms represent ensembles, like molecules or clusters.

Example:

    UnitCell[
      GdFeVoid_site=Variant [
        (p=1/3)
        Gd  1   0 0 0        0.035 0.035 0.035 0 0 0
        
        (p=1/3)
        iron_molecule = [
          Fe 1 0 0  0.2845   0.035 0.035 0.035 0 0 0
          Fe 1 0 0 -0.2845   0.035 0.035 0.035 0 0 0
        ]
        
        (p=1/3)
        Void 
      ]
    ]

In this example `GdFeVoid_site` is either occupied by a gadolinium atom, by a structural building block of two iron atoms, which is aliased as `iron_molecule`, or by a `Void`. On average, each chemical unit is present with the same probability of `p=1/3`.

### Groups of atoms

Atoms can be grouped into one entity. This entities represent atoms which can not appear without each other, usually molecules or clusters.

The groups are defined without keyword, just with square brackets `[` `]`.

Example:

    Fe2_molecule = [
      Fe 1   0 0  0.2845   0 0 0 0 0 0
      Fe 1   0 0 -0.2845   0 0 0 0 0 0
    ]

The grouping is recurrent, each group may contain not only atoms, but also other groups.

Example:

    Fe2_molecule = [
      
      Fe 1   0 0  0.2845   0 0 0 0 0 0

      lower_iron = [ 
        Fe 1   0 0 -0.2845   0 0 0 0 0 0
      ]
    ]


### Atoms

Atoms can be defined in two formats with both isotropic and anisotropic ADP:

    AtomType  mult  x y z  Uiso
    AtomType  mult  x y z  U11 U22 U33 U12 U13 U23

`x` `y` and `z` are fractional coordinates, `Uiso` and `Uij` are atomic displacement parameters in Å$^2$. The multiplier `mult` can be thought of as a site multiplicity in the average structure. The occupancy of each atoms is defined as `mult*p` where `p` is the probability of the atomic group the atom belongs to (defined in the construction `(p=...)`) and `mult` is the multiplier. In the vast majority of the cases `mult` should be equal to 1. It is only different if the atomic group contains disorder which is not disentangled into different entries of `Variant`.

Example:

    Fe   1   0 0 -0.2845   0 0 0 0 0 0

### Naming the entities of the unit cell

Variants, atomic groups and atoms can be given names for symmetry expansion and use in `Correlations` and `Modes`. Example:

    Variant1 = Variant[
    	(p=1)
    	Iron_molecule = [
    		Fe1 = Fe 1   0 0  0.2845   adp_perp adp_perp adp_z 0 0 0
    		Fe2 = Fe 1   0 0 -0.2845   adp_perp adp_perp adp_z 0 0 0
    	]
    ]

### Symmetry

It is possible to apply symmetry elements to atoms and groups of atoms. The syntax is the following:

     atomic_group*Symmetry(x,y,z)

Example for mirroring an atom with respect to the xy0 plane:

     Fe_atom = Fe 1   0 0  0.2845   0.01
     
     symmetry_equivalent_Fe = Fe_atom*Symmetry(x,y,-z)


Unlike in the average crystal structure, the symmetry elements are sensitive to integer translations. In cases when the center of the molecule is not in `(0,0,0)` user should provide the appropriate translations. Here is example for the same iron molecule, with centered at `(0,0,1)`

    Fe_atom = Fe 1   0 0  0.7155   adp_perp adp_perp adp_z 0 0 0
    
    symmetry_equivalent_Fe = Fe_atom*Symmetry(x,y,2-z)    #note 2-z here
    # resulting coordinates are 0 0 1.2845



## Modes 

Modes is a way to define motions of atoms and molecules in Yell. All the motions are linear meaning that a displacement of each atom $\Delta r$ is expressed as a linear combination of translation along all the modes:

$$
\Delta r_i = \sum_n M^n_i \xi_n
$$

where $M^n_i$ is a set of vectors attached to each atom in a molecule, $\xi_n$ is the amplitude of a displacement of the molecule along mode $n$.

Currently Yell supports two types of modes: translations and linear approximation of rotations.

### TranslationalMode
Defines movements of atomic groups as a whole; expressed in Å.

Format: `TranslationalMode(atomic_group,axis)` 
or `TranslationalMode(atom,axis)`
where `atomic_group` is a name to a group of atoms, and `axis` could be either `x`,`y` or `z`. 

Example:
    
    manganese_x = TranslationalMode(manganese,x)

### RotationalMode

Defines linear approximation to rotation of a group of atoms along arbitrary axis; expressed in radians.

Format: `RotationalMode(atoimc_group, axis_x, axis_y, axis_z, center_x,center_y,center_z)`  
where `axis_x` `axis_y` and `axis_z` define the rotation axis direction and `center_x` `center_y` `center_z` define a point on the rotation axis. 

A displacement vector for each atom is defined as 

$$
r_{axis} \times \frac{r_{atom}-r_{center}}{|r_{atom}-r_{center}|}.
$$ 

## Correlations

In this section, short-range order pairwise correlations are defined.  There are three basic types of correlations: substitutional correlations, displacement correlations and size effect (see [Weber & Simonov, 2012](#WeberSimonov)). The combination of these correlations may describe any disorder in a disordered crystal.

In Yell, correlations which belong to the same group of atoms should be combined together in groups using square brackets. Example:

    Correlations [
      [(1,0,0)
       Multiplicity 2
       SubstitutionalCorrelation(CuVoid,AuVoid,0.25+dp)
       SizeEffect(Cu,Au_x,dx)
       ADPCorrelation(Cu_x,Au_x,dUx)
       ]
     ]
     
Each group starts with a vector `(u,v,w)`, defining the lattice vector between the correlated atoms or molecules. It corresponds to the $\mathbf R_{u,v,w}$ in ([Weber & Simonov, 2012](#WeberSimonov)). The elements of the vector are not required to be integer, and could for example be multiples of 0.5 (in case of C centering) and even irrational (e.g. in case of quasicrystals). The corresponding inter-atomic vector is calculated as $\mathbf R_{u,v,w} + \mathbf r_{Au} - \mathbf r_{Cu}$ , where $\mathbf r_{Cu} $ and $\mathbf r_{Au}$ are the atomic coordinates of the copper and gold atoms defined within `CuVoid` and `AuVoid` variants (not shown in the example).

### Multiplicity {#Multiplicity}
Format: `Multiplicity m`  
where `m` is number.

Defines the multiplicity of the interatomic or intermolecular pair with respect to the Laue group.

Effect: multiplies `m` with the multiplicity of each atomic pair in current correlation group. Does not have a letter dedicated to it in ([Weber & Simonov, 2012](#WeberSimonov)), but is equivalent in effect by multiplication `m` to both joint probability $p^{mn}_{uvw}$ and Patterson multiplicity $c_nc_m$.

In Yell, the Laue group is applied to the diffuse scattering automatically. However, Yell currently cannot automatically calculate the multiplicity of each correlation set. It is expected that the user manually provides the multiplicity for it.

The multiplicity `m` of a group should be equal to the number of times such correlation group appears in the ∆PDF space. The correlations in the center of ∆PDF space usually get multiplicity 1, the correlations on general positions get multiplicity that is equal to the order of the crystal Laue group.

For an example see the model for **FeVoid**.

**Warning:** The definition of the multiplicity has to be done with great care. Wrongly defined multiplicities may be perfectly compensated by over- or underestimated structure correlation parameters leading to errors that are difficult to recognize!

For detailed instructions on how to apply multiplicity see [this section](#HowMultiplicity).


### Substitutional correlation

Substitutional correlation means that the occupancies of two `Variants` in the structure depend on each other.

Effect: sets the joint probabilities of atomic pairs to the ones specified as arguments. Affects $p^{mn}_{uvw}$ as stated in ([Weber & Simonov, 2012](#WeberSimonov)).

Format: `SubstitutionalCorrelation(variant_A,variant_B,p11,p12,...,pnm)` 
short format: `SubstitutionalCorrelation(variant_A,variant_B,p11,p12,...,pn_minus_1m_minus_1)`  
where `pij` are the joint probability coefficients to find both blocks Ai and Bj present at the same time, the `pn_minus_1m_minus_1` is the joint probability $p_{n-1,m-1}$.

#### Example:

    SubstitutionalCorrelation(CuAu,CuAu,0.25+x)

#### Extended example

For the unit cell:

    UnitCell [
      AuCu = Variant[
        (p=1/3)
        Au 1  0 0 0  0
        (p=2/3)
        Cu 1  0 0 0  0
    ]]
    Correlations [
      [(0,0,0)
      SubstitutionalCorrelation(AuCu,AuCu,1/3) #corresponds to the matrix 1/3 0
                                               #                          0 2/3
      ]
    ]

This will produce the following pairs (differences are typed in red):

<pre>
      m p        x y z  Uxx ... Uyz  p̅        x̅ y̅ z̅  U̅xx ... U̅yz
Au Au 1 <font color="brown">0.333333</font> 0 0 0  0 0 0 0 0 0  0.111111 0 0 0  0 0 0 0 0 0
Cu Au 1 <font color="brown">0</font>        0 0 0  0 0 0 0 0 0  0.222222 0 0 0  0 0 0 0 0 0
Au Cu 1 <font color="brown">0</font>        0 0 0  0 0 0 0 0 0  0.222222 0 0 0  0 0 0 0 0 0
Cu Cu 1 <font color="brown">0.666666</font> 0 0 0  0 0 0 0 0 0  0.444444 0 0 0  0 0 0 0 0 0
</pre>

Here the `m` marks pair multiplicicty, `p` is joint probability, `x` `y` `z` are interatomic vector coordinates, `Uxx ... Uyz` are the components of the joint ADP tensors, the letters with overbar `p̅ x̅ y̅ z̅  U̅xx ... U̅yz`  relate to the average probability, interatomic vector and ADP tensor.

#### Formal description

Assume that 

    variant_A = Variant [
      (p=pA1)
      A1
      (p=pA2)
      A2
      ...
      (p=pAn)
      An
    ]
    variant_B = Variant [
      (p=pB1)
      B1
      (p=pB2)
      B2
      ...
      (p=pBm)
      Bm
    ]
Where `Ai` and `Bj` are some chemical units, or void. Then, the full matrix of joint probabilities has size $m \times n$  and the form

$$
P = \left[
\begin{array}{cccc}
p_{A1B1} & p_{A2B1} & ... & p_{AnB1} \\
p_{A1B2} & p_{A2B2} & ... & p_{AnB2} \\
... & ... && ... \\
p_{A1Bm} & p_{A2Bm} & ... & p_{AnBm} \end{array} \right]
$$

This will be expressed in yell in the following way:

    SubstitutionalCorrelation(variant_A,variant_B,pA1B1,pA2B1,...,pAnB1,
                                                  pA1B2,pA2B2,...,pAnB2,
                                                  ...
                                                  pA1Bm,pA2Bm,...,pAnBm)
                                                  
Since possibilities `Ai` are all exclusive and $\sum_i p_{Ai}=1$ there are $m+n-1$ independent constraints on $p_{AiBj}$ which constrain the sum of pair probabilities to probabilities of single chemical units:

$$
\sum_i p_{AiBj} = p_{Bj}
$$
$$
\sum_j p_{AiBj} = p_{Ai}
$$

The constraints make it possible to calculate the last row and last column of the matrix $P$, resulting in

$$
P = \left[
\begin{array}{cccc}
p_{A1B1}  & ... & p_{An-1B1} & p_{B1}- \sum_{i=1}^{n-1} p_{AiB1}  \\
... & ... && ... \\
p_{A1Bm-1}  & ... & p_{An-1Bm-1} & p_{Bm-1} - \sum_{i=1}^{n-1} p_{AiBm-1} \\
p_{A1} - \sum_{j=1}^{m-1} p_{A1Bj} & ... 
& p_{An-1} - \sum_{j=1}^{m-1} p_{An-1Bj} & p_{AnBm}
\end{array} \right]
$$
where
$$
p_{AnBm}=1+\sum_{i,j=1}^{n-1,m-1}p_{Ai,Bj}-\sum p_{Ai}-\sum p_{Bj}  
$$

Yell can automatically calculate the last row and the last column, given the upper left $(m-1)\times(n-1)$ independent part of joint probability matrix. Thus the above definition can be equally well expressed in the short form:

    SubstitutionalCorrelation(variant_A,variant_B,pA1B1,pA2B1,...,pAn-1B1,
                                                  pA1B2,pA2B2,...,pAn-1B2,
                                                  ...
                                                  pA1Bm-1,pA2Bm-1,...,pAnBm-1)
The short form is in general preferred, because it is a non-redundant representation of the structure model.

#### Neutral joint probabilities

If the occupancies of two variants are independent, e.g. because their interatomic distance is larger than the correlation length of the corresponding local order phenomenon, the joint probabilities are equal to the product of the occupancies $p_{AiBj}^{Patterson}= p_{Ai}p_{Bj}$. Thus the following probability matrix will produce no diffuse scattering and no signals in PDF space (and can be ommited):

$$
P_{neutral} = \left[
\begin{array}{cccc}
p_{A1}p_{B1} & p_{A2}p_{B1} & ... & p_{An}p_{B1} \\
p_{A1}p_{B2} & p_{A2}p_{B2} & ... & p_{An}p_{B2} \\
... & ... && ... \\
p_{A1}p_{Bm} & p_{A2}p_{Bm} & ... & p_{An}p_{Bm} \end{array} \right]
$$

#### Zero - neighbor joint probabilities

For the zero neighbor correlation, i.e. the correlation of a `Variant`with itself, the diagonal elements of the joint probability matrix should be equal to the average occupancy of the corresponding elements, and off-diagonal elements should be equal to zero:

$$
P_{zero} = \left[
\begin{array}{cccc}
p_{A1} & 0 & ... & 0 \\
0 & p_{A2} & ... & 0 \\
... & ... && ... \\
0 & 0 & ... & p_{An} \end{array} \right]
$$

Or in the Yell (long) form:

    [(0,0,0)
     SubstitutionalCorrelation(variant_A,variant_A,pA1,0,...,0,
                                                   0,pA2,...,0,
                                                   ...
                                                   0,0,...,pAn)]
It is important to note that the zero-neighbor correlation **must** be defined for all Variants showing occupational disorder and that the joint probabilities should not be refined but be the same as the occupation probability (p=...) defined in the UnitCell statement.  

In practice the zero neighbor correlation is frequently over- or underestimated, because incomplete background subtraction or elimination of broad diffuse scattering with the background (= overcorrection of background) may influence the zero neighbor correlation definition and thus the correct determination of the scale factor. As shown in [(Weber & Simonov, 2012)](#WeberSimonov), such errors may lead to significant systematic deviations in the determination of other pair correlation parameters. Careful background determination is therefore indispensable for a high quality local structure refinement.

## ADPCorrelation
This correlation appears when displacements of two atoms or molecules are not independent.

Effect: changes $\beta^{mn}_{uvw}$ as stated in ([Weber & Simonov, 2012](#WeberSimonov))

Format: `ADPCorrelation(Mode1,Mode2,cov)`
where `cov` is a covariance of the displacements along the two modes $cov = < \xi_1 \xi_2 >$; `cov` is expressed in $\text{units}_1*\text{units}_2$ where $\text{units}_i$ are the units of the two modes (e.g. Å for translational modes and rad for rotational ones).

The correlations of atomic displacements typically manifest themselves as thermal diffuse scattering (TDS), but the correlations could also be of static origin. 

Yell assumes that displacements of all atoms in the crystal are jointly Gaussian and the distribution of the displacement differencess $\Delta r_{ik}-\Delta r_{jl}$ is also Gaussian. The correlations of the displacements are expressed in terms of covariances of collective displacement of the blocks along `Modes`. Analogously to independent substitutional correlations ADPCorrelations do not need to be defined, if `Modes` are not correlated with any other mode. In fact, most of the internal modes of rigid molecules are not expected to be correlated.

Note that some ADP correlations, though symmetrically independent, can produce the same signal in PDF space. One very practical example is the covariances $<x_1y_2>$ and $<y_1x_2>$ produce exactly the same effect in PDF space, and should therefore be constrained. 

The ADP correlations are subject of symmetry constraints by the symmetry of interatomic or intermolecular pairs. The constraints to the correlations of translational modes are equivalent to the site-symmetry ADP constraints in the average structure. The constraints involving rotational modes are equivalent to the TLS constraints ([Schomaker & Trueblood 1968](#ShomakerTrueblood)).

### Example

    #Cu structure (fcc)
    #the unit cell has Fm-3m symmetry (225)
    LaueSymmetry m-3m
    UnitCell [
      Variant [
        (p=1)
        Cu = Cu 1  0 0 0  0.01
    ]]
    Modes [
      Cu_x = TranslationalMode(Cu,x)
    ]
    Correlations[
    ...
    [(1,0,0)
     ADPCorrelation(Cu_x,Cu_x,0.001)
     ADPCorrelation(Cu_y,Cu_y,0.0001)
     ADPCorrelation(Cu_z,Cu_z,0.0001)]
    ...
    ]

Note that the yy and zz correlations are the same, while xx is different according to symmetry of the pair. Let us assume that the symmetry of this pair is described by the point group $4/mmm$ with the 4-fold axis along **a** of the cubic crystal. This restricts the covarince matrix to the following form:

$$
\left[
\begin{array}{ccc}
cov_{xx} &  &   \\
 & cov_{yy} &   \\
 &  & cov_{yy}
\end{array} 
\right]
$$

### Example with correlation matrices

The covariance matrices can be recalculated in correlation matrices. 

$$
corr^{mn}_{ij} = \frac{cov^{mn}_{ij}}{\sqrt{\mathbf U^{aver,m}_{ii} \mathbf U^{aver,n}_{jj}}}
$$

where indices $i,j=1,2,3$ mark tensor components, $m$ and $n$ count atoms,  $corr^{mn}_{ij}$ is a correlation matrix, $cov^{mn}_{ij}$ is a covariance matirx (in Å$^2$), the $U^{aver,m}$ and $U^{aver,n}$ are the average ADP tensors of atoms $m$ and $n$ (in Å$^2$).

In Yell correlation matrix can be calculated in the epilogue using the `Print` command. The above mentioned example becomes:

    #Cu structure (fcc)
    #the unit cell has Fm-3m symmetry (225)
    LaueSymmetry m-3m
    UnitCell [
    Uiso=0.01; #define a variable for Uiso
      Variant [
        (p=1)
        Cu = Cu 1  0 0 0  U
    ]]
    Modes [
      Cu_x = TranslationalMode(Cu,x)
    ]
    Correlations[
    ...
    cov_11=0.001;
    cov_22=0.002;
    [(1,0,0)
     ADPCorrelation(Cu_x,Cu_x,cov_11)
     ADPCorrelation(Cu_y,Cu_y,cov_22)
     ADPCorrelation(Cu_z,Cu_z,cov_22)]
    ...
    ]
    Print "corr_11=" cov_11/Uiso " corr_22=" cov_22/Uiso " corr_33=" cov_22/Uiso

Will print:

    Requested output: corr_11=0.1 corr_22=0.01 corr_3=0.01




### Formal definition {#ADPCorr}

The theory is somewhat similar to the normal mode analysis. The notation here is analogous to ([Cyvin 1968](#Cyvin)).

 **1.** Denote the equilibrium configuration of N atoms in a molecule $\alpha$ as
 
 $$
 \boldsymbol{R}_{\alpha i}=\{X_{\alpha i}, Y_{\alpha i}, Z_{\alpha i} \}=\{R_{\alpha i1}, R_{\alpha i2}, R_{\alpha i3} \}
 $$

for i=1,2,...,N. Then the displacements are introduced as deviations from equilibrium:

$$
u_{\alpha ij}=R_{\alpha ij}-<R_{\alpha ij}>
$$

**2.** Introduce another set of coordinates $\boldsymbol \xi_{\alpha}$ which we identify as a set of internal coordinates and introduce a set of matrices $M$ that transform internal coordinates into "external" atomic displacements:

$$
u_{\alpha ij}=\sum_k M^k _{\alpha ij} \xi_{\alpha k}
$$

If $\xi$ has 3N independent coordinates and there exist an inverse transformation

$$
\xi_{\alpha i}=\sum_{jk} M^{-1 jk} _{\alpha i} u_{\alpha jk}
$$

both representations are equivalent.

In Yell, the matrices $M^k _{\alpha ij}$ are called *Modes*, the variables $\xi_{\alpha i}$ are referred to as *the amplitudes of displacement along corresponding modes*.

**3.** Recall, that the elements of the ADP matrix of an individual atom $k$ in the molecule $\alpha$ is defined as follows:

$$
U_{\alpha k}^{ij} = <u_{\alpha k i} u_{\alpha k j}>
$$

The paper ([Weber & Simonov, 2012](#WeberSimonov)) uses the notation $\beta_{\alpha k}^{ij}$; it is equivalent notation since $\beta_{\alpha k}^{ij}=2\pi^2a^*_i a^*_j U_{\alpha k}^{ij}$, where $a^*_i$ and $a^*_j$ are the lengths of the corresponding reciprocal lattice vectors.

**4.** A covariance between the displacements of two atoms belonging to two different molecules can be expressed in terms of covariances of the modes:

$$
<u_{\alpha i j} u_{\beta k l}>=<\sum_m M^m_{\alpha i j}\xi_{\alpha m} \sum_nM^n_{\beta k l}\xi_{\beta n}>=\sum_{mn} M^m_{\alpha i j} M^n_{\beta k l} <\xi_{\alpha m} \xi_{\beta n}>
$$

**5.** The joint atomic displacement parameter matrix of one atom as seen from another atom can be expressed as follows:

$$
U_{\alpha k,\beta l}^{ij} = <(u_{\alpha k i}-u_{\beta l i})(u_{\alpha k j}-u_{\beta l j})>=\\
<u_{\alpha k i}u_{\alpha k j} >+<u_{\beta l i} u_{\beta l j}>-
<u_{\beta l i}u_{\alpha k j}>-<u_{\alpha k i}u_{\beta l j}>=\\
U_{\alpha k}^{ij}+U_{\beta l}^{ij}-\sum_{mn} (M^m_{\alpha i j} M^n_{\beta k l}+M^n_{\alpha i j} M^m_{\beta k l}) <\xi_{\alpha m} \xi_{\beta n}>
$$

Note that in order to define the ADP correlations it is not necessary to construct the full set of internal modes. It is significant to know the ADPs of each pf the atoms from the average structure and only the covariances of the modes which are correlated.

### Neutral correlation

If the displacements of two molecules are independent, all the covariances are equal to zero and don't need to be listed explicitly.

### Zero-neighbor correlation

In theory, the covariance in the zero neighbor is equal to the square of the amplitude of the mode. For example, for single atom, the ADP correlations of the zero neighbor should be identical with its ADP parameter: 


    adp_cu=0.01;
    
    UnitCell [
      Variant [
        (p=1)
        Cu = Cu 1  0 0 0  adp_cu
    ]]
    Modes [
      Cu_x = TranslationalMode(Cu,x)
      Cu_y = TranslationalMode(Cu,y)
      Cu_z = TranslationalMode(Cu,z)
    ]
    Correlations [
      [(0,0,0)
       ADPCorrelation(Cu_x,Cu_x,adp_cu)
       ADPCorrelation(Cu_y,Cu_y,adp_cu)
       ADPCorrelation(Cu_z,Cu_z,adp_cu)]
    ]

Similar to substitutional correlations zero neighbor ADP correlations may be heavily affected by over- or under-corrected experimental background.

Note that in substitutionally disordered crystals ADP correlations are frequently not seen because corresponding diffuse scattering is usually very weak. In such crystals zero-neighbor ADP correlation need only be defined if there is a clear indication to ADP correlations or size-effect correlations.

## SizeEffect

Size effect occurs when a substitutional disorder in a `Variant` is correlated with a static displacement of another chemical unit.

Effect: changes the $u^{mn}_{uvw}$  as stated in ([Weber & Simonov, 2012](#WeberSimonov))
Format: `SizeEffect(AtomicGroup1,Mode2,amp)` 
or `SizeEffect(Mode1,AtomicGroup2,amp)` 

The meaning of the `SizeEffect(AtomicGroup1,Mode2,amp)`: in the presence of `AtomicGroup1` the `AtomicGroup2` has an average displacement along the `Mode2` by the amount defined by the amplitude `amp` (expressed in the units of the `Mode2` e.g. Å for translational mode, rad for rotational)

### Example

    UnitCell [
      Variant [
        (p=0.99)
        Na = Na  1  0 0 0  0.01
        
        (p=0.01)
        Void
    ]]
    Modes [
      Na_x = TranslationalMode(Na,x)
    ]
    Correlations [
      ...
      [(1,0,0)
       SizeEffect(Na,Na_x,0.02)
      ]
      ...
    ]

### Neutral SizeEffect

When the correlation is absent, the size effect is equal to zero and is not required to be defined. 

### Zero-neighbor correlation

The size-effect for the zero-neighbor is equal to zero. However, it is important to add the zero neighbor from the ADPCorrelation, otherwise the model might produce negative intensities.

Example:

    #A disordered Fe-Ni alloy
    #Space group Fm-3m
    PointGroup m-3m
    UnitCell
    [
    ADP=0.0004;
      FeNi = Variant[
        (p=1/2)
        Fe = Fe 1  0 0 0  ADP
    
        (p=1/2)
        Ni = Ni 1  0 0 0  ADP
      ]
    ]
    
    Modes[
      Fe_x = TranslationalMode(Fe,x)
      Ni_x = TranslationalMode(Ni,x)
      #... same for y and z
    ]
    Correlations [
      [(0,0,0)
       SubstitutionalCorrelation(FeNi,FeNi,1/2)
       #IMPORTANT: even though the crystal does not have ADP correlation,
       #zero neighbor ADP correlation should be input in presense of size effect:
       ADPCorrelation(Fe_x,Fe_x,ADP) 
       ADPCorrelation(Ni_x,Ni_x,ADP)
       #... same for y and z
      ]
    #...
    #
      [(1,0,0)
       SubstitutionalCorrelation(FeNi,FeNi,0.25*(1+a200))
       
       SizeEffect(Fe,Fe_x,se200_FeFe)
       #IMPORTANT: in addition to size effect add ADP correlation 
       #with amplitude se^2/2
       ADPCorrelation(Fe_x,Fe_x,pow(se200_FeFe,2)/2) 
       
       SizeEffect(Fe,Ni_x,se200_FeNi)  
       ADPCorrelation(Fe_x,Ni_x,pow(se200_FeNi,2)/2)
       
       SizeEffect(Ni,Fe_x,se200_FeNi)
       ADPCorrelation(Ni_x,Fe_x,pow(se200_FeNi,2)/2)
       SizeEffect(Ni,Ni_x,se200_NiNi)
       ADPCorrelation(Ni_x,Ni_x,pow(se200_NiNi,2)/2)
      ]
    

### Formal definition

In the presence of size-effect:

$$
\bar u^{\alpha k,\beta l}_i = \sum_nM^n_{\beta l i}A_{\beta n}-\sum_nM^n_{\alpha k i}A_{\alpha n}
$$

For description of the notation see [ADPCorrelation](#ADPCorr)

Epilogue
------------

In the last part of the file the custom output from Yell can be requested. The output will be produced after the refinement is finished. The output can contain any expressions whose values will be calculated from the refined variables.

### Print
optional.
Arguments: a list of strings in quotes `"` `"` or arithmetic expressions
Format: `Print "string1" expr1 ...`

Prints the requested output in the terminal.

Example:
    
    a=1;b=2;
    Print "The variable a=" a ", b=" b 

This line produces the output `Requested output:The variable a=1, b=2`

How to calculate multiplicity {#HowMultiplicity}
=========

Yell can apply the Laue symmetry to the calculated ∆PDF. This is a convenient feature since only the [independent cone](#Cone) of correlations in ∆PDF space have to be defined. 

Yell does not distinguish the interatomic pairs on special and average positions. Thus, the multiplicity of each interatomic (intermolecular) pair **must be provided manually**. If multiplicity is not provided, all the pairs except in the center of ∆PDF will get too little density:

![enter image description here][1]

If the multiplicity is provided, the Laue symmetry is applied properly:

![enter image description here][2]

Definition
---------

By *multiplicity* of an interatomic (intermolecular) pair we mean the number of interatomic (intermolecular) pairs which are symmetry related to the current pair. 

By symmetry we mean two distinct symmetries: the crystal space group symmetry and [combinatorial symmetry](#Combinatorial) which relates pair (A,B) to (B,A).

First way to determine multiplicity
-----------

Multiplicity can be determined by counting the symmetry equivalent pairs. 

### Example 
This example is performed in two dimensions, but generalization to three dimensions is straightforward.

Assume a simple two-dimensional square crystal with a space group $p4mm$ and one atom in the unit cell:

![][3]

The zeroth neighbor connects an atom with itself. There is 1 such pair:

![][4]

The first neighbor with interatomic vector (1,0) has 4 symmetry equivalents:

![][5]

The neighbor (1,1) also has 4 symmetry equivalents:

![][6]

also the neighbor (2,0):

![enter image description here][7]

The neighbor (2,1) has 8 equivalents:

![enter image description here][8]

The neighbor (2,2) has again 4 equivalents:

![enter image description here][9]

The list of all pairs (without actual correlations) in Yell format will be the following:

    Correlations [
      [(0,0,0)
       Multiplicity 1
       ...]
       
       [(1,0,0)
       Multiplicity 4
       ...]
       
       [(1,1,0)
       Multiplicity 4
       ...]
       
       [(2,0,0)
       Multiplicity 4
       ...]
       
       [(2,1,0)
       Multiplicity 8
       ...]
       
       [(2,2,0)
       Multiplicity 4
       ...]
     ]
       
Second way to determine multiplicity {#FromSymmetry}
-----------

Sometimes it is complicated to count the number of symmetry related pairs or there is a demand to calculate the multiplicity by an alternative approach to verfy the results. In such cases, the multiplicity can be calculated from the internal symmetry of a single pair.

Assume:

* $p$ is a pair
* $|G|$ is the number of symmetry elements in the crystal space group (space group order)
* $|Int(p)|$ is the number of symmetry elements in the internal symmetry of a pair (internal symmetry order)

Then, number of symmetry equivalent pairs is equal to:

$N_p = 2\frac{|G|}{|Int(p)|}$ if the pair connects different atoms( molecules) and
$N_p = \frac{|G|}{|Int(p)|}$ for zeroth neighbor.

The above mentioned formula is a consequence of the [orbit-stabilizer theorem](http://www.proofwiki.org/wiki/Orbit-Stabilizer_Theorem). The factor 2 appears in the formula due to [compinatorial symmetry](#Combinatorial). 

### Example

Again, assume a simple square crystal. The crystal has a plane group $p4mm$ (No. 11 in International Tables of Crystallography):

![enter image description here][10]

The plane group $p4mm$ has following symmetry operations:

$$
\begin{array}{llll}
\text{(1) $x,y$} & \text{(2) $\bar x,\bar y$} & \text{(3) $\bar y,x$} & \text{(4) $y,\bar x$} \\
\text{(5) $\bar x,y$} & \text{(6) $x,\bar y$} & \text{(7) $y,x$} & \text{(8) $\bar y,\bar x$}
\end{array}
$$

Totally there are 8 operations, thus $|G|=8$.

The zeroth neighbor has internal symmetry $4mm$:

![enter image description here][11]

The order of internal symmetry group $4mm$ is 8 (group elements: 1, $4^1$, $4^2$, $4^3$, $m_x$, $m_y$, $m_{xy}$, $m_{yx}$); this pair is a zeroth neighbor:

$$
N_{(0,0)} = \frac{|G|}{|Int(p_{(0,0)})|} = \frac{8}{8}=1
$$

Hence there is 1 such pair.

The neighbor (1,0) has internal symmetry $2mm$:

![enter image description here][12]

On the image the glide planes which are marked orange belong to the space group of the crystal, but not to the internal symmetry of the pair.

The order of $2mm$ is 4 (group elements: 1, 2, $m_x$, $m_y$), thus the number of such pairs is:

$$
N_{(1,0)} = 2\frac{|G|}{|Int(p_{(1,0)})|} = 2\frac{8}{4}=4
$$

Note that glide plane operations do no count in this case, because the internal symmetry of the pair has to be described exclusively by point symmetry operations. In the case of periodic layers or rods, which are disordered along the other dimensions (corresponding to diffuse streaks or layers, respectively), glide or screw operations might also be considered along the periodic directions (see corresponding layer or rod groups). 

The neighbor (1,1) also has internal symmetry $2mm$

![enter image description here][13]

The average crystal has two additional planes and a four fold axis (marked orange) passing through coordinates (0.5,0.5), which do not completely belong to the internal symmetry of the pair. They are not counted. The remaining group elements are : 1, $4^2\equiv 2$, $m_{xy}$, $m_{yx}$

The number of (1,1) pairs is therefore:

$$
N_{(1,1)} = 2\frac{|G|}{|Int(p_{(1,1)})|} = 2\frac{8}{4}=4
$$

Neighbor (2,0) again has the symmetry $2mm$

![enter image description here][14]

and $N_{(2,0)} = 2*8/4=4$

Neighbor (2,1) has the symmetry $2$ (group elements: 1, 2):

![enter image description here][15]

Thus the number of such pairs is $N_{(2,0)} = 2*8/2=8$

Neighbor (2,2) has symmetry $2mm$ (group elements: 1, 2, $m_{xy}$, $m_{yx}$):

![enter image description here][16]

And thus there are 4 such pairs.

From symmetry analysis it is clearly seen:

* The zeroth neighbor has symmetry $4mm$ and only 1 equivalent.
* Neighbors $(x,0)$ and $(x,x)$ have symmetry $2mm$ and thus 4 equivalents.
* All the other neighbors $(x,y)$ have symmetry $2$ and thus 8 equivalents.

Combinatorial symmetry {#Combinatorial}
-----------

By *combinatorial equivalent* of an (ordered) pair (A,B) we mean the pair (B,A). In PDF, combinatorial pairs are obtained by inverting interatomic vectors.

Assume we have a planar NaCl-type structure:

![enter image description here][17]

The pair Na-Cl with interatomic vector (0.5,0.5) has a combinatorially equivalent pair Cl-Na with vector (-0.5,-0.5). Thus the multiplicity of such pair is 8:

![enter image description here][18]



When the multiplicity is calculated using  [the second method](#FromSymmetry), combinatorially symmetric pairs are automatically counted. The pair Na-Cl (0.5,0.5) has the internal symmetry $m$:

![][19]

The multiplicity is equal to $N_{NaCl} = 2*8/2=8$

Independent cones for all Laue groups {#Cone}
-----------

The following table summarizes the independent parts of each Laue group:

$$
\begin{array}{ccc}
\text{Laue group} &\text{Group order} &\text{Independent cone conditions} \\
m\bar 3m& 48& x \geq y \geq z \geq 0 \\
m\bar 3 & 24 & x \geq z,\; y \geq z,\; z \geq 0 \\
6/mmm & 24 & x \geq 2y \geq 0,\; z\geq0 \\
6/m & 12 & x\geq y\geq 0,\; z \geq 0 \\
\bar 3m:H &  12  &  x \geq y \geq 0,\; z \geq 0 \\
\bar 3m:R &  12  &  z \geq y \geq x,\; x+y+z \geq 0 \\
\bar 3:H &  6  & x \geq 0,\; y \geq 0,\; z \geq 0 \\
\bar 3:R & 6  & x \geq y,\; x \geq z ,\; x+y+z \geq 0 \\
4/mmm  &  16  & z \geq 0,\; x \geq y \geq 0 \\
4/m  &  8  &  x \geq 0,\; y \geq 0,\; z \geq 0 \\
mmm   &  8  &  x \geq 0,\; y \geq 0,\; z \geq 0 \\
2/m  &  4  &  z \geq 0,\; y \geq 0 \\
2/m:b  &  4  & z \geq 0,\; y \geq 0 \\
\bar 1  &  2  &  z \geq 0
\end{array}
$$


Equation for diffuse scattering calculation {#TheFormula}
=============

Diffuse scattering of a disordered crystal can be calculated by the equation (8) from ([Weber & Simonov 2012](#WeberSimonov):

$$
I_{dif}(\mathbf{h}) = 
\sum_{\mathbf{R}_{uvw}}^{cryst} \sum_{mn}^{cell} 
\big\{p^{uvw}_{mn}
\exp(-\mathbf{h}^{T}\beta^{mn}_{uvw} \mathbf{h})
\cos  [2\pi \mathbf h (\mathbf{R}_{uvw}+\mathbf{r}_{mn}+ \mathbf{\bar u}^{mn}_{uvw} ]-\\
- c_{m}c_{n} \exp(-\mathbf{h}^{T}(\beta^{aver}_m+\beta^{aver}_n) \mathbf{h}) \cos  [2 \pi \mathbf{h}(\mathbf{R}_{uvw} + \mathbf{r}_{mn}) ] 
\big \}  \cdot \\ 
 \cdot f_{m}(\mathbf{h}) f_{n}(\mathbf{h}) 
$$

Here indices $m$ and $n$ go over all the atoms in a unit cell, $\mathbf{R}_{uvw}$ over all latice vectors; $p^{uvw}_{mn}$ is a joint probability to find atom $n$ in a unit cell and atom $m$ in another unit cell separated by $\mathbf{R}_{uvw}$, $\beta^{mn}_{uvw}$ is a joint ADP matrix expressed in fractional units, $r_{mn}=r_n-r_m$ is the average interatomic vector between atoms $m$ and $n$, $\mathbf{\bar u}^{mn}_{uvw}$ is a size effect parameter, $c_m$ and $c_n$ are the average occupancies, and $\beta^{aver}_m$ and $\beta^{aver}_n$ are the average ADP parameters, $f_m(\mathbf h)$ and $f_n(\mathbf h)$ are the atomic form-factors.

References
========

##### Weber, T., & Simonov, A. (2012). The three-dimensional pair distribution function analysis of disordered single crystals: basic concepts. *Zeitschrift für Kristallographie*, 227(5), 238-247. [researchgate](https://www.researchgate.net/publication/235767166_The_three-dimensional_pair_distribution_function_analysis_of_disordered_single_crystals_basic_concepts)  {#WeberSimonov}


##### Schomaker, V., & Trueblood, K. N. (1968). On the rigid-body motion of molecules in crystals. *Acta Crystallographica Section B: Structural Crystallography and Crystal Chemistry*, 24(1), 63-76. [link](http://scripts.iucr.org/cgi-bin/paper?S0567740868001718)  {#ShomakerTrueblood}


##### Cyvin, S. J. (1968). *Molecular vibrations and mean square amplitudes.* Universitets Forl.. Chapter 3  {#Cyvin}


##### Simonov, A., Weber, T., Steurer, W., *Yell - a computer program for diffuse scattering analysis via 3D-∆PDF refinement.* in. prep.  {#YellPaper}



  [1]: https://lh6.googleusercontent.com/ijtDa2ou82uPzYmpkhAoJwkRg-uYx2teR2AKzBVyzu0=s600
  [2]: https://lh6.googleusercontent.com/1iW1RNDe3IXG6yjvmxF44UPsS_RfR28LoXJi5Do7BYc=s600
  [3]: https://lh3.googleusercontent.com/8YG3ZmWPOcUOr5B6iMYc8aE408ETOHXlng2q1Xv_4zM=s300
  [4]: https://lh5.googleusercontent.com/029euAKAuzzT0_EjpQR36ghC7Z1kfWfxM4SIBq9St7o=s300
  [5]: https://lh4.googleusercontent.com/mHDYn0VHH8Q9VIaNy3mBbZXb8cv0oPHb-YStyhumpl8=s300
  [6]: https://lh5.googleusercontent.com/GH8IkDiOtM4DfR5scDmrvQwRjhz7OM8kcmZf00mUcaQ=s300
  [7]: https://lh3.googleusercontent.com/Guw5EK2cBAsKsPH4_MNabz5Plg5vL4s3erYF5Y11FB8=s300
  [8]: https://lh6.googleusercontent.com/CBZwrTbE3CtRohPa9u1L7Y4Gr-Qlq3jrP4VpbsU9GSY=s300
  [9]: https://lh3.googleusercontent.com/XZrzH28SR3XrQSl2nJlB0aoUkflqMM9bV1XfL93Dlo4=s300
  [10]: https://lh4.googleusercontent.com/_q2pPCcpWjywdCnL391OhM7GeeiT7X9nnLCRsGIYBhA=s600
  [11]: https://lh3.googleusercontent.com/QH2LQbeljajEbxPlJlc5lmaa7bssClaQ--4pO3gxuyA=s600
  [12]: https://lh6.googleusercontent.com/p2Y_0uI0QJ-E4sfOTb_oMpHPghPcogAU7H7DlMM3RiM=s600
  [13]: https://lh5.googleusercontent.com/XPNTxskn1FCCARzcPrYk6JZbs696oMoWr0VeiixeZ40=s600
  [14]: https://lh4.googleusercontent.com/Xx0l08giXHN-azKwZQlFE7joCcg0U0-AGHsKXehGCBE=s600
  [15]: https://lh5.googleusercontent.com/xYBJYzKWc1meOpN1KcmHq77KihMfuMGMuSdEufRNFFw=s600
  [16]: https://lh3.googleusercontent.com/wHBFjQQuDMyJhj0ZizRvT0fvtpndqoTIIQBdZyAnarc=s600
  [17]: https://lh5.googleusercontent.com/AIBvkOR4j_t9K1sr0HifYynBOGsd59xy29vopSox2h0=s300
  [18]: https://lh4.googleusercontent.com/YcF870Bckp7yCAvJPswlwTb1jstT6OkFaCwTOWl9fH4=s300
  [19]: https://lh4.googleusercontent.com/-GHrkdfqmj8jFQbvINVpzuTfppmu8Md6aR3uhoiRxvI=s450