import numpy as np
import itertools
import os
import shutil

"""Script for preparing a large scanning table and splitting it to smaller subscans.
"""

def makeRangeInclusive(start: float, stop: float, stepSize: float) -> np.ndarray:
    """Uses smaller step size in the last step if the division is not even"""


    if stop <= start:
        return np.array([start])

    # This is sensitive to floating point errors, the division can be eg. 5.99999 => we get 5, which is wrong. We "fix" this below
    wholeIntervals = int((stop - start) / stepSize)

    values = [start + i*stepSize for i in range(wholeIntervals)]

    ## Now append enough full step sizes to ensure we reach stop. If we go past it, undo one step and append stop directly 
    while values[-1] < stop:
        values.append( values[-1] + stepSize )

    if values[-1] > stop:
        values[-1] = stop

    return np.array(values)


def main():

    """
    values_mh2 = makeRangeInclusive(100, 800, 5)
    values_a2 = makeRangeInclusive(0, 8, 0.1)
    values_b4 = makeRangeInclusive(0.25, 2.0, 0.25)
    values_sintheta = makeRangeInclusive(-0.3, 0.3, 0.02)
    values_b3 = makeRangeInclusive(-150, 150, 50)
    """
    values_mh2 = makeRangeInclusive(100, 150, 50)
    values_a2 = makeRangeInclusive(0, 1, 0.2)
    values_b4 = makeRangeInclusive(0.25, 0.25, 0.25)
    values_sintheta = makeRangeInclusive(-0.1, 0.1, 0.1)
    values_b3 = makeRangeInclusive(-50, 50, 50)
    
    ## ---- Construct all possible combinations of these and put them in a 2D table

    scanTable = np.array(list( itertools.product(values_mh2, values_a2, values_b4, values_sintheta, values_b3) ))

    numPoints = len(scanTable)

    tableHeader = ["mh2", "a2", "b4", "sinTheta", "b3"]

    ## ---- Prepare N directories, one for each subscan, and copy necessary files there
    numProcesses = 1

    print(f"Preparing {numProcesses} subscans, {int(numPoints / numProcesses)} points each. {numPoints} points total.")

    ## Should we copy and prepare files for slurm batch submission?
    bSlurmFiles = True

    filesToCopy = ["parameters"]

    if bSlurmFiles:
        filesToCopy.append("submit.job")
        # We give each job a unique name, so read the file here for easier replacing later
        with open("submit.job", "r") as f:
            jobFileLines = f.readlines()
        
        for i, line in enumerate(jobFileLines):
            if line.startswith("#SBATCH --job-name"):
                jobNameLineIndex = i


    def copyFiles(dirName: str) -> None:
        # helper function
        for fname in filesToCopy:
            pathToCopy = shutil.copy(fname, dirName)

            if fname == "submit.job":
                # TODO should not actually need to copy this since we have the contents stored
                jobFileLines[jobNameLineIndex] = (f'#SBATCH --job-name="scan{i+1}"')
                with open(pathToCopy, "w") as f:
                    f.writelines(jobFileLines)


    ## ---- Slice the parameters table into a smaller 2D table and store in appropriate subdir
    numRows, numColumns = scanTable.shape
    rowsPerDir = numRows // numProcesses
    remainder = numRows % numProcesses

    for i in range(numProcesses):
        dirName = f"subscan_{i+1}"
        os.makedirs(dirName, exist_ok=False)

        copyFiles(dirName)

        startRow = i * rowsPerDir
        endRow = startRow + rowsPerDir + (1 if i < remainder else 0)

        scanTableSlice = scanTable[startRow:endRow, :]

        filePath = os.path.join(dirName, "scanningTable.csv")
        np.savetxt(filePath, scanTableSlice,  delimiter=",", header=",".join(tableHeader), fmt='%g')


if __name__ == "__main__":
    main()