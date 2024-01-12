# Miniroot

Read root file storing TTree (single but can be easily extended to multiple), assuming the types of all branches are trivial and known.

## Problems

- How to retrieve type information of each branch (from TTree Metadata at the end? Try to save the uncompressed TTree metadata)

> Difficulty of deserializing TTree Metadata, may change between versions
