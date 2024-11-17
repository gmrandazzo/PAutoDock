#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""__autogridmap2dx__.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

Provides a commandline script to recover the output result of a virtualscreening.
"""

import argparse
import sys


class AutoGridMap2DX:
    def __init__(self, map_file=None):
        self.name = ""
        self.npts = [0, 0, 0]
        self.n = [0, 0, 0]
        self.center = [0, 0, 0]
        self.origin = [0, 0, 0]
        self.nelem = 0
        self.spacing = 0.0
        self.values = []
        self.datafile = ""
        self.molecule = ""
        self.paramfile = ""
        self.precision = 0.0001
        if map_file is not None:
            self.name = map_file
            with open(map_file, "r") as fp:
                self.read(fp)

    def read(self, fp):
        for i in range(6):
            line = fp.readline()
            if i == 0:
                self.paramfile = line.split()[1]
            elif i == 1:
                self.datafile = line.split()[1]
            elif i == 2:
                self.molecule = line.split()[1]
            elif i == 3:
                self.spacing = float(line.split()[1])
            elif i == 4:
                self.npts = [int(x) for x in line.split()[1:]]
            elif i == 5:
                self.center = [float(x) for x in line.split()[1:]]
        for i in range(3):
            self.n[i] = self.npts[i] + 1
        self.nelem = self.n[0] * self.n[1] * self.n[2]
        i = 0
        while i < self.nelem:
            val = float(fp.readline())
            self.values.append(val)
            i += 1
        for i in range(3):
            self.origin[i] = self.center[i] - self.npts[i] / 2 * self.spacing

    def writeDX(self, fname):
        with open(fname, "w") as fp:
            nx, ny, nz = self.n
            ori = self.origin
            spacing = self.spacing
            vals = self.values

            print("#==================================", file=fp)
            print(f"# AutoGrid Map File: {self.name}", file=fp)
            print(f"# Receptor File Name: {self.molecule}", file=fp)
            print("#==================================", file=fp)
            print(f"object 1 class gridpositions counts {nx} {ny} {nz}", file=fp)
            print(f"origin {ori[0]:12.5E} {ori[1]:12.5E} {ori[2]:12.5E}", file=fp)
            print(f"delta {spacing:12.5E} {0:12.5E} {0:12.5E}", file=fp)
            print(f"delta {0:12.5E} {spacing:12.5E} {0:12.5E}", file=fp)
            print(f"delta {0:12.5E} {0:12.5E} {spacing:12.5E}", file=fp)
            print(f"object 2 class gridconnections counts {nx} {ny} {nz}", file=fp)
            print(
                f"object 3 class array type double rank 0 items {len(vals)} data follows",
                file=fp,
            )
            col = 0
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        fp.write(f" {vals[i*ny*nz + j*nz + k]:12.5E}")
                        col += 1
                        if col == 3:
                            print(file=fp)
                            col = 0
            if col != 0:
                print(file=fp)
            print('attribute "dep" string "positions"', file=fp)
            print('object "regular positions regular connections" class field', file=fp)
            print('component "positions" value 1', file=fp)
            print('component "connections" value 2', file=fp)
            print('component "data" value 3', file=fp)


def main():
    """
    main.py
    """
    p = argparse.ArgumentParser()
    p.add_argument("--map", default=None, type=str, help="autodoc4 map")
    p.add_argument("--dx", default=None, type=str, help="DX output")
    args = p.parse_args(sys.argv[1:])

    if args.map is None or args.dx is None:
        print("\nUsage: %s --map [input autodock4 map]" % sys.argv[0])
        print("                --dx [output DX]")
    else:
        agm = AutoGridMap2DX(args.map)
        agm.writeDX(args.dx)


if __name__ in "__main__":
    main()
