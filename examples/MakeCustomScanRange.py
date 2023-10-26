#!/Users/Osmanthus/opt/anaconda3/bin/python3
import numpy as np
import argparse

## This script prepares an example input file for SingletEWPT 
# that can be used for non-uniform parameter scanning (logarithmic example included)

## Prepare range of values so that their 10-base logarithms are evenly spaced
def makeLogRange(min: float, max: float, npoints: int) -> np.ndarray:

    try: 
        start = np.log10(min)
        end = np.log10(max)
        exponentRange = np.linspace(start, end, npoints)
        return 10**exponentRange
    except Exception as e: 
        print(e)
    
## Just a uniform grid
def makeUniformRange(min: float, max: float, npoints: int) -> np.ndarray:
    return np.linspace(min, max, npoints)

## Random 
def makeRandomRange(min: float, max: float, npoints: int) -> np.ndarray:
    randomRange = np.random.uniform(min, max, npoints)
    return np.sort(randomRange)

def makeInputFile(paramName: str, values: np.ndarray) -> None:
    fileName = "range_" + paramName
    with open(fileName, 'w') as file:
        for x in values:
            file.write(f'{x}\n')


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('name', type=str, metavar="Parameter name", choices=['mh2', 'sinTheta', 'a2', 'b3', 'b4', 'T'])
    parser.add_argument('min', type=float, metavar="Starting value (inclusive)")
    parser.add_argument('max', type=float, metavar="Ending value (inclusive)")
    parser.add_argument('npoints', type=int, metavar="Number of points")
    parser.add_argument('--type', type=str, metavar="Type of value range to use. Uniform by default.", default='uniform', choices=['uniform', 'log10', 'random'])

    args = parser.parse_args()

    if (args.type == 'uniform'):
        values = makeUniformRange(args.min, args.max, args.npoints)
    elif (args.type == 'log10'):
        values = makeLogRange(args.min, args.max, args.npoints)
    elif (args.type == 'random'):
        values = makeRandomRange(args.min, args.max, args.npoints)

    makeInputFile(args.name, values)


if __name__ == "__main__":
    main()