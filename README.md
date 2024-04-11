<!-- <a name="readme-top"></a> -->
<!-- Table of Contents -->
<details>
<summary>Table of Contents</summary>

- [CASToR Reconstruction for HiRezBrainPET](#castor-reconstruction-for-hirezbrainpet)
  - [Content](#content)
    - [*BPET\_castor\_v3.1.1* folder](#bpet_castor_v311-folder)
    - [*BPET\_castor\_results* folder](#bpet_castor_results-folder)
  - [Getting Started](#getting-started)
    - [Pre-requisites](#pre-requisites)
    - [Installation](#installation)
  - [How to use](#how-to-use)
    - [CASToR Source code](#castor-source-code)
    - [Results - Python scripts](#results---python-scripts)
  - [Roadmap](#roadmap)
  - [License](#license)

</details>

<!-- END Table of Contents -->

<!-- CASToR Reconstruction for HiRezBrainPET -->
# CASToR Reconstruction for HiRezBrainPET

[CASToR](https://castor-project.org/) is an open source software for tomographic reconstruction, developed by the LaTIM laboratory in France, which stands for Customizable and Advanced Software for Tomographic Reconstruction. With flexibility, modularity and generic programming in mind, the main aim is creating a frame work for developing and testing new reconstruction methods and algorithms.

This repository contains modifications to the code and new the scripts to run the reconstruction program for the HiRezBrainPET system, an innovative PET scanner that uses RPC as detectors. The [first results](https://arxiv.org/abs/2211.05860) were published in November 2022 by Paulo Fonte et al. using a simple reconstruction process. Thus, the goal with this is to make that process more generic and customizable to the system, and allowing a more flexible and scalable tool for researchers in the field of medical imaging, or even the broader public, to use and replicate the results. All of this in a more efficient and easier way for the development, implementation and testing.

## Content

The project consists of 2 main folders, the main source code of the CASToR project with the modifications needed to account for the HiRezBrainPET system, ``BPET_castor_v3.1.1``, and python program, the batch scripts and the graphical user interface (using tkinter from python) to run the program for the HiRezBrainPET system, `BPET_castor_results`.

### *BPET_castor_v3.1.1* folder

The main changes lie in the `toolkits` folder, which consist of the CASToR utilities that helps of the data processing. The ones relevant and modified/created are:

- `castor-PETScannerLutEx_BPET.cc` - The main program to create Look-Up-Table (LUT) files (.lut) that represents any scanner geometry, like the HiRezBrainPET system
- `castor-scannerLUTexplorer.cc` - A program to visualize the LUT file content, exploring the central positions and orientation of the approximated system crystal, element by element.
- `castor-txtConversionCrystalsID_BPET.cc` - A program to convert the LOR coordinates to the crystal indexation in the scanner geometry. It returns a similar txtfile with the format t, id1, id2, where t is the time of detection, and id1 and id2 are the crystal id of both intersection points, instead of the x, y, z coordinates.
- `castor-datafileConversionEx_BPET.cc` - A program to convert the txt file to the CASToR format.

Its with these programs that the data is prepared to be used in the main reconstruction program, `castor-recon.cc`.
After one needs to compile the project, following the instructions in section [Usage](#usage). MODIFIED WITH CORRECT SECTION

### *BPET_castor_results* folder

Here we have a **configuration** folder, which can be deprecated if the one used in the build folder, BPET_castor_v3.1.1, is used.

Additionally, we have a **Release** folder with the the executables of the programs, including the main reconstruction program, `castor-recon.exe`, as well as the batch files to run the programs with the correct parameters. There are several batch files to run the main program, all with the same purpose, named as `runmpi_2024Derenzo_win{}.bat`, where the {} is the number of a certain test for specific parameters.
In addition, there is a subfolder **batch_python** with a code for a graphical user interface (GUI) using tkinter, that creates the batch files given the user input, that can be saved and run the reconstruction program with the corresponding parameters.

Finally, there is a **python** folder with the scripts to process the data of the results, that adapt the data to the format of the reconstruction program, and to visualize the results in a more intuitive way. The instructions on how to use the scripts are in the [Usage](#usage) section. MODIFIED WITH CORRECT SECTION

## Getting Started

### Pre-requisites

For the python scripts, it is compatible with all operating systems that support Python 3.11 or above.

As for the compilation of the CASToR source code, it is necessary to have the following installed:

- CMake 3.21 or newer
- .NET Core 8.0 or newer - for Windows
- MSVC 2015 (Visual Studio) or newer - for Windows
- GCC 4.8 or newer - for Unix-based systems
- If using MPI (Message Passing Interface) for parallel processing, it is necessary to have the MPI library installed, both the standalone executables and the SDK installers.

An example of the paths needed to be added/check in the environment variables, in Windows:

- `C:\Program Files\CMake\bin` - CMake
- `C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\Bin` - MSVC
- `C:\Program Files\Microsoft MPI\Bin` - MPI
- `C:\Program Files (x86)\Microsoft SDKs\MPI\Include\` - MPI
- `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\` - MPI
- `C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\` - MPI
- `C:\Program Files\dotnet\` - .NET Core

### Installation

Starting with the python scripts, it is necessary to install the requirements in the `requirements.txt` file, found in the python subfolder of the BPET_castor_results folder. This can be done with the following command:

```bash
pip install -r requirements.txt
```

As for the CASToR source code, a guide can be found in section 3 of the CASToR Documentation, [here](https://castor-project.org/documentation_v3). The main steps though are:

Create a build folder in the root of the project and run cmake to build the project, using the following commands:

```bash
mkdir build
cd build
cmake ..
```

Alternatively, one can use the **CMake GUI**, `cmake-gui` to configure the build, selecting the source code folder and the build folder, and then clicking on the `Configure` button. After that, one can click on the `Generate` button to create the build files using the either the MSVC or the GCC compiler or even other compiler depending on the operating system.

Once created the makefile for the compiler, one can run the following command to compile the project, in case of using the GCC compiler:

```bash
make
```

If using the MSVC compiler, one can open the solution file created in the build folder, and build the project using the Visual Studio IDE. Make sure to build the project in Release mode, to have the best performance.
However, there is a way to compile the project directly in the command line, using the msbuild command, as follow:

```bash
msbuild CASTOR.sln /p:Configuration=Release /t:ALL_BUILD
```

where the `/p:Configuration=Release` flag is to build the project in Release mode, and the `/t:ALL_BUILD` flag is to build all the targets in the project. Other target projects can be available to build, which can be seen in the build folder.

## How to use

### CASToR Source code

Here there are several changes that can be done to the source code to adapt to the HiRezBrainPET system, depending on what is needed.

The main changes currently made are in the `toolkits` subfolder, as mention in section [Content](#content), with the creation/adaptation of the programs: `castor-PETScannerLutEx_BPET.cc`, `castor-txtConversionCrystalsID_BPET.cc` and `castor-datafileConversionEx_BPET.cc`.

Beyond the creation and adaption of tools to process the data, CASToR also allows the creation of new modules, by building specific classes that inherits from more abstract classes. In CASToR [website](https://castor-project.org/documentation_v3), there are guides on how to create new modules, as well as specific documentation for specific type of classes, like:

- image convolver (such as PSF model)
- image processing
- projectors
- optimizers and penalizers
- dynamic models, etc...

All of this module added in both source (src) and header (include) folders.

As for the main program, `castor-recon.cc`, it can also be adapted, although it is not very recommended, as it is the main program that runs the reconstruction process. However, it is possible to add new parameters to the program, or even new options to the program, by adding new flags to the program, and then adapting the code to read and use these flags.

> **IMPORTANT:** After each change in the source code, it is necessary to recompile the project, as mentioned in the [Installation](#installation) section. For example, if a toolkit is changed/created, it is necessary to recompile the project to have the new executable of that program (In case of the creation of a new program, the CMakeList.txt file must be updated to include the new program in the build process).

### Results - Python scripts



In order to run the program open a terminal window and run the following command:

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

config folder of the results should be the same as the one in the build folder.

## Roadmap

## License

[MIT](https://choosealicense.com/licenses/mit/)
