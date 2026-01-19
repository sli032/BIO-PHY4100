# viral capsid modeling

# Version
gsd 2.9.0
hoomd 3.5.0

# Compile hoomd-blue
Compile hoomd following the steps:
1. Download and extract [**hoomd 3.5.0**](https://github.com/glotzerlab/hoomd-blue/releases/tag/v3.5.0).
2. Download 'bondflip_plugin' folders and put them inside 'hoomd-v3.5.0/example_plugins/'.
3. Open 'hoomd-v3.5.0/example_plugins/CMakeLists.txt' and add 'add_subdirectory(bondflip_plugin)' at the end of file.
4. Open terminal and type
```bash
$conda create -n hoomd35 python=3.10
$conda activate hoomd35
$conda install cmake eigen=3.4.0 pybind11

$cd hoomd-blue-3.5.0
$mkdir build
$cd build
$cmake ..
$make -j20
$make install
```
*Note for new mac user, replace files "hoomd/extern/nano-signal-slot/nano_signal_slot.hpp" and "hoomd/GPUFlags.h".\

5. Now compile and install plugins:
```bash
$cd ..
$cmake -B build/example_plugins -S example_plugins
$cmake --build build/example_plugins
$cmake --install build/example_plugins
```

*Folder 'bondflip_plugin' is developed for hoomd-blue 3.5.0. 
For a general plugin template and use, please refer to [Plugins in HOOMD-Blue](https://hoomd-blue.readthedocs.io/en/v3.5.0/components.html).


# Python script
The python script for virus assembly is included in the pythonexample folder.
Before you run the python script, download the necessary packages:
```bash
conda install numpy pandas pyyaml
conda install -c conda-forge gsd=2.9.0
```

## Set up parameter
The parameter space is included in "parentparam.yml", the parameter could be either ON/OFF, or string or a list of values.

## Generate folder with given parameters
Run 
$ py param.py parentparam.yml.

It will iterate all parameters and generate folders with json file included.

Copy files inside 'messyT=3' folder to the new generated parameter folder. (in this way you will start from a messy geometry with the same number of triangles as a T=3 capsid)

## Run simulation
$./runscript.sh "folder name"

If permission is required, run

$ chmod +x runscript.sh

before run the .sh file.

# Visualization
To visualize the simulation, download [**Ovito**](https://www.ovito.org), and open the 'color.gsd' file inside the result folder.
