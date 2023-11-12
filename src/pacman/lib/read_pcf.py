"""This class loads a PACMAN control file (pcf) and lets you
   querry the parameters and values.

   Constructor Parameters
   ----------------------
   file : pathlib.Path
       A control file containing the parameters and values.

   Notes
   -----
   A parameter can have one or more values, differet parameters can
   have different number of values.

   The function Param.get(index) automatically interprets the type of the
   values. If they can be cast into a numeric value retuns a numeric
   value, otherwise returns a string.

   Examples
   --------
   # Load a pcf file
   >>> import reader3 as rd
   >>> reload(rd)
   >>> pcf = rd.Pcffile('/home/patricio/ast/esp01/anal/wa011bs11/run/wa011bs11.pcf')

   # Each parameter has the attribute value, wich is a ndarray:
   >>> pcf.planet.value
   ...array(['wa011b'], dtype='|S6')

   # To get the n-th value of a parameter use pcffile.param.get(n):
   # if it can't be converted to a number/bool/etc, it returns a string.
   >>> pcf.planet.get(0)
   ... 'wa011b'

   >>> pcf.photchan.get(0)
   ... 1

   >>> pcf.fluxunits.get(0)
   ... True

   # Use pcffile.param.value[n] to get the n-th value as string:
   >>> pcf.aorname.get(0)
   ... 38807808

   >>> pcf.aorname.value[0]
   ... '38807808'

   # The function pcffile.param.getarr() returns the numeric/bool/etc
   # values of a parameter as a nparray:
   >>> pcf.sigma.value
   ... array(['4.0', '4.0'], dtype='|S5')

   >>> pcf.sigma.getarr()
   ... array([4.0, 4.0], dtype=object)

   Modification History
   --------------------
   2009-01-02 chris Initial Version.
                    by Christopher Campo      ccampo@gmail.com
   2010-03-08 patricio Modified from ccampo version.
                       by Patricio Cubillos      pcubillos@fulbrightmail.org
   2010-10-27 patricio Docstring updated
   2011-02-12 patricio Merged with ccampo's tepclass.py
   2021-12    Sebastian Zieba Updated for PACMAN usage
"""
from pathlib import Path
from typing import Any

import numpy as np

from .options import OPTIONS


# NOTE: Each parameter is an instance of this class
class Param:
    def __init__(self, vals: Any) -> None:
        """The class's constructor."""
        self.value = vals

    def get(self, index: int = 0) -> Any:
        """Return a numeric/boolean/None/etc. value if possible,
        else return a string."""
        try:
            return eval(self.value[index])
        except Exception:
            return self.value[index]

    def getarr(self):
        length = np.size(self.value)
        ret = np.zeros(length, dtype='object')
        for i in np.arange(length):
            ret[i] = self.get(i)
        return ret


class Pcf:
    def __init__(self, params: np.ndarray) -> None:
        """The class's constructor."""
        for parname in params:
            value = None
            if "dir" in parname[0]:
                value = [Path(parname[1]).resolve()]
            else:
                value = parname[1:]
            setattr(self, parname[0], Param(value))

    def make_file(self, name: Path) -> None:
        with name.open('w', encoding=OPTIONS["encoding"]) as file:
            attrib = vars(self)
            keys = attrib.keys()

            file.write(f"@ {self.pcfname.get()}\n")
            for key in keys:
                if key != "pcfname":
                    file.write(f"{key} {attrib.get(key).value[0]}\n")


def read_pcf(file: Path) -> None:
    """Function to read the file."""
    # NOTE: List containing the set of parameters:
    pcfsets = []

    # NOTE: Read the file
    with file.open('r', encoding=OPTIONS["encoding"]) as f:
        lines = f.readlines()

    cleanlines = []     # List with only the important lines
    block = []          # Blocks separator
    # NOTE: Clean the lines:
    for i in np.arange(len(lines)):
        line = lines[i]

        # NOTE: Strip off comments:
        try:
            line = line[0:line.index('#')].strip()
        except Exception:
            line = line.strip()

        # NOTE: Keep only useful lines:
        if len(line) > 0:
            cleanlines.append(line)
            # NOTE: Identify the separators:
            if line[0] == "@":
                block.append([len(cleanlines) - 1, line[1:].strip()])

    # NOTE: Append a line to mark the end of the block:
    block.append([len(cleanlines), "end"])

    if (len(block) - 1) == 0:
        # NOTE: Do normal readpcfs
        params = []
        for line in cleanlines:
            params.append(np.array(line.split()))
        return Pcf(params)
    else:
        # NOTE: Loop over each block:
        for i in np.arange(len(block) - 1):
            params = []     # List for the parameters and values of the block
            multiples = []  # Line position of multiple valued parameter
            nval = []       # Number of values

            for j in np.arange(block[i][0] + 1, block[i + 1][0]):
                params.append(np.array(cleanlines[j].split()))

                # NOTE: If the parameter has more than 1 value:
                if len(params[-1]) > 2:
                    multiples.append(len(params) - 1)
                    nval.append(len(params[-1]) - 1)

            # NOTE: Number of parameters with multiple values:
            nmul = len(multiples)

            if nmul == 0:
                pcfsets.append(params)
                pcfsets[-1].append(["pcfname", str(block[i][1])])
            else:
                # NOTE: Calculate total number of sets
                nt = 1
                for j in np.arange(nmul):
                    nt *= nval[j]
                ncurrent = nt

                # NOTE: Holder of the sets of params
                # and make nt copies of the original set
                parset = []
                for j in np.arange(nt):
                    parset.append(params[:])

                    # NOTE: Add the pcfname:
                    parset[j].append(["pcfname", str(block[i][1])])

                # NOTE: Loop over each multiple valued parameter:
                for j in np.arange(nmul):
                    ncurrent /= nval[j]
                    mpar = np.copy(params[multiples[j]][1:])

                    # NOTE: Edit the value in each set:
                    for k in np.arange(nt):
                        index = int((k / ncurrent) % nval[j])
                        parset[k][multiples[j]] = np.array(
                                [params[multiples[j]][0], mpar[index]])
                for ps in parset:
                    pcfsets.append(ps)

        # NOTE: Return a List of Pcf objects (one for each set):
        return [Pcf(pcfset) for pcfset in pcfsets]


def store_pcf(meta, pcf: Pcf) -> None:
    """Store values from PACMAN control file as parameters in Meta object."""
    for key in pcf.__dict__.keys():
        try:
            setattr(meta, key, pcf.__dict__[key].get(0))
        except Exception:
            try:
                setattr(meta, key, pcf.__dict__[key].getarr(0))
            except Exception:
                print("Unable to store parameter: ", key)
