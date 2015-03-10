SETTING UP:

do the first two steps if the libfftw3-3.lib file does not already exist in the
"fftw3" directory

1. 	Open up a [visual studio 2010] command line terminal.

2. 	change directory ("cd") into the "ext" project folder and enter ther following command: 
	"C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\amd64\lib.exe" /def:libfftw3-3.def /MACHINE:X86
_____________________________________________________________________________
3. 	create a new directory "build" (or whatever you what it to be)

4. 	"cd" into that directory and invoke the cmake command as follows: 
	cmake -G "Visual Studio 10" ..

5. 	open your windows explorer and locate your current directory within the terminal. In there should be
	the visual studio solution file. Open it.

6. give me an "A" 

	B-)