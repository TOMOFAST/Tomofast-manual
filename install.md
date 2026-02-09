## Tomofast-x Install

Installation instructions for Windows and Linux systems are provided below.

- [Windows Native](#windows-native)
- [Windows WSL](#windows-wsl)
- [Linux / MacOS](#linux--macos)


#### Windows Native (runtime, Recommended!)
1) Download and install oneAPI runtimes:   
   
- https://tectonique.net/tomofast-x-q/intel-mpi-2021.17.2.93_offline.exe   and 
     
- https://tectonique.net/tomofast-x-q/w_ifx_runtime_p_2025.3.2.835.exe      


2) Download the precompiled tomofast-x executable:
   
- https://tectonique.net/tomofast-x-q/tomofastx.exe     

3) Open a **Command Prompt** console from **Start Menu** (click on start icon then type **cmd** and the Command Prompt tool will be shown)

- Once the console is open, type (including quotes where shown):
``` 
“C:\Program Files (x86)\Intel\oneAPI\setvars.bat”
``` 
4) Then change directory (linux command is cd) to downloaded tomofastx.exe directory, then type:
``` 
mpiexec -n 4 tomofastx.exe -p parfiles\Parfile_mansf_slice.txt
```  


#### Windows Native (for Compilation)
1) Install Visual Studio Build Tools and C++ Desktop Tools   
   
- https://visualstudio.microsoft.com/downloads/?q=build+tools   
   
- Select to install BOTH Desktop development with C++ AND C++ Tools for Linux and Mac Development   
    
2) Install IntelOne API Toolkit   

- https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html

- Follow Instructions

3) Optionally Install GNU Make for Windows [Optional, if you intend to compile tomofast-x rather than using a precompiled binary]   

- https://gnuwin32.sourceforge.net/packages/make.htm

- Follow Instructions



##### To compile the source code from scratch:
- i.	Download and unzip latest tomofast-x code from: 

- https://github.com/TOMOFAST/Tomofast-x

- ii.	Open **x64 Native Tools** Command Prompt from Start Menu → Visual Studio 2026 Directory

- iii.	Once the console is open, copy or type (including quotes):

``` 
 “C:\Program Files (x86)\Intel\oneAPI\setvars.bat”
``` 
- iv.	Then copy or type:

``` 
set PATH=%PATH%;C:\Program Files (x86)\GnuWin32\bin
``` 
- v.	Then change directory (linux command is cd so something like cd C:\Users\vogarko\Downloads\Tomofast-x-master) to unzipped Tomofast-x code directory, then copy or type:
``` 
make WINDOWS=1
``` 
- vi.	If the code compiles without error, it will create a new tomofastx.exe file and then you can test the code with:
``` 
mpiexec -n 4 tomofastx.exe -p parfiles\Parfile_mansf_slice.txt
```

#### Windows WSL
First install Windows Subsytem for Linux (WSL): Open **Windows PowerShell** or Windows **Command Prompt** in administrator mode by right-clicking and selecting **Run as administrator**, type in the following command on the command line  
``` 
    wsl --install   
```
then restart your machine. This will also install the Ubuntu operating system wihin WSL (it will not affect your normal Windows system). Now follow the steps for Linux installation below.

#### Linux / MacOS
1) Install **gfortran** compiler using the relevant install commands: https://fortran-lang.org/learn/os_setup/install_gfortran/   
2) Install **OpenMPI** library using the relevant install commands: https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html   
3) If not already installed, install **git** on your computer using the relevant install commands: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git   
4) Download **Tomofast-x** code by entering: 
```
    git clone https://github.com/TOMOFAST/Tomofast-x.git  
```

5) Change to the tomofast-x directory and compile the code by entering the following command:   
```
     make
```
6) If the code compiles without error, it will create a new tomofastx file and then you can test the code with:   
```
     mpirun -np 4 ./tomofastx -p ./parfiles/Parfile_mansf_slice.txt
```

Note that the Makefile by default assumes the gfortran (GCC) compiler. It can be switched to the Intel compiler by setting the flag "COMPILER = 2". To use the Intel compiler, you will need to install it separately on your system.
