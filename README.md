This is the source codes for the Calibur program reported in https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-25

(For precompiled Windows binary of Calibur, please visit https://sourceforge.net/projects/calibur/)

How to compile:

   If you are running a variant of UNIX (or MacOS/WSL/Cygwin/MSYS2), you should have the GNU C compiler (`gcc`) and Make (`make`) installed.
   At the directory where the file Makefile is, type `make`.
   This will build two programs: calibur and calibur-lite
   
   - calibur is the main program.
   - calibur-lite is a lighter version of calibur which reports only the best decoy.

   If you are running MS Windows, you should have installed MS Visual Studio Community.
   Open the Developer Command Prompt at the directory where the file Makefile.win is, type `nmake /F Makefile.win`.
   This will build the program calibur.exe.

   That's all!

How to run:

1. Prepare a file which list all the decoys' PDB files.
   That is, the file should read like:
```
decoysdir/1enhA0.pdb
decoysdir/1enhA1000.pdb
decoysdir/1enhA1001.pdb
decoysdir/1enhA1002.pdb
decoysdir/1enhA1003.pdb
...
```
   Alternatively, you can use Rosetta's `compose_score_silent.py` to create a
   silent file out of your PDB files, and use this file as input to calibur.

2. Call calibur with this filename (`mydecoys` say) as the only parameter

   `sh> calibur mydecoys`

   calibur will figure out a reasonable threshold itself. You can also specify
   a threshold. For example, for threshold 1.0, run calibur as

   `sh> calibur mydecoys 1.0`

   You can also change the strategy that calibur uses to discover the
   threshold. To see how, run

   `sh> calibur`

   And the options will be shown.
