# fastmerge
The fastmerge tool has been developed to quickly merge together molecules saved as .js files made in VMD.
The usual use case for this would be if multiple seperate files have been prepared for size reasons, but need to be combined at the file level.
I have also found it useful when reordering atoms found in polymer fragments.

Basic compilation instructions are to follow the makefile, which compiles the relatively simple executable which the user can move into their path.
```bash
make
mv fastmerge /somewhere/in/the/path
```
To use the tool, the command can take any number of arguments, with the first being the output file and the rest being inputs.

```bash
fastmerge outputfile.js input1.js input2.js ...
```
