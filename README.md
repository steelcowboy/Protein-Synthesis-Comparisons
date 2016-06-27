# Protein Synthesis Comparisons v0.4
Programs to convert DNA to RNA and find protein sequences,
implemented in C++ and Python to compare speeds. 

The C++ version is completely usable with many different
flags. The Python version will eventually be updated to 
match the functionality of the C++ version.

## To Use
For general purposes use cmake. To keep everything clean:

```bash
mkdir _build
cd _build
cmake ../
make
```

To test the programs against each other 
`cd src/comparison` run `make`.

*Note:* Make comparison will create the protein_synthesis 
binary without any cout statements

## To Do
- ~~Move everything to a separate library~~
- Allow input from a file
- Integrate Python into the build system
- See how my Python and C can go together (Cython?)
- Clean up some implementation
  * Separate output into sections and allow querying certain information
- Bring Python version up-to-date
- Make documentation



