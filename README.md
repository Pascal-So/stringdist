# stringdist - Map out strings on a 2d plane

This program takes a list of strings and returns a list of 2d coordinates, one coordinate for every string in the input, representing a possible embedding of the levenshtein distances between the strings in the 2d plane, using sammon mapping.

NOTE: This code is slooooow.. building the distance matrix is O(n^2 * l^2), where n is the number of strings and l is the length of the longest string, but it's fast enough to have some fun embedding strings in to the plane :)

## Usage

prepare a file with one string per line, and pass it in to the standart input of stringdist. The first and only argument to stringdist is the amount of passes for the gradient descent, default is 500 if no value is passed.

Example:
```bash
cat names.txt | ./stringdist 1000 | gnuplot -p -e "p '-'"
```
