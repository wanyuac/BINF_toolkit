"""
Convert MUMmer's output into a CSV file.
Yu Wan (20 Apr 2017)
Example: python tabulateMUMerCoordinates.py input.coords > output.coords
Reference: David Edwards' script filterCoords.py in the RedDog (https://github.com/katholt/RedDog) suite.
Licence: Apache-2.0
Python version: 3.5.2 (but compatible to Python 2)
"""

import sys

def main():
    count = 0
    with open(sys.argv[1], "rU") as f:
        print("Start,End,Identity")
        for line in f:
            if count <= 5:  # skip the first five lines, including the self-self comparison (100% identity)
                count += 1
            else:
                data = line.split("|")
                identity = float(data[3])  # drops all white spaces
                coords = data[0].split()  # removes all white spaces as well
                start = coords[0]
                end = coords[1]
                print(",".join([start, end, str(identity)]))

if __name__ == "__main__":
    main()