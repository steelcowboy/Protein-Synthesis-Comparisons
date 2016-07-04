# 0.5
- Made a few other optimizations, mostly regarding adding some const methods
- Fixed bug where one sequence without a start sequence crashes the program

# 0.4
- Moved functions into a library
- Optimized a few things
- Added new class members to avoid passing arguments:
  * organism_number
  * start_sequence
  * end_sequence
  * is_sequenced
- Added new function to retrieve sequenced value
- Added the number argument to the Organism constructor
- Moved all output related operations to the get_output function

# 0.3
- Switched to getopt for argument parsing
- Added a default option which has an example DNA set
