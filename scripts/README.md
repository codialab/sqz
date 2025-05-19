# Scripts

## compress-gfa-sequitur.sh
This script can be used to compress a GFA file using [Sequitur](https://github.com/craignm/sequitur) (we used the C++ implementation). For it to work the program `sequitur` has to be available on the path.

To run use:
```bash
./scripts/compress-gfa-sequitur.sh <GFA-FILE> 2
```
(`2` means the number of appearances when a new `Q`-line should be created)

## count-path-length.py
Used to get the lengths of paths in nodes for the analysis of `sqz`.

## coverage.py
Can be used to calculate the coverage of nodes in either a compressed or uncompressed GFA file. Used for the analysis of `sqz`
