#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mgltoolsinstall.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

Provides the basic operation for molecular files.

"""

import logging
import platform
import subprocess
from pathlib import Path

from pautodock.fileutils import download_file, extract_tar_gz


def install_mgltools(dir: str):
    """
    Install MGLTools from the official website.
    """
    directory = Path(dir)
    directory.mkdir(exist_ok=True)
    destination = directory / "mgltools_1.5.7.tar.gz"
    dir_name = None

    system = platform.system()
    if system == "Linux":
        url = "https://ccsb.scripps.edu/download/532/"
        dir_name = "mgltools_x86_64Linux2_1.5.7"
    elif system == "Darwin":
        url = "https://ccsb.scripps.edu/download/529/"
        dir_name = "mgltools_1.5.7_MacOS-X"
    else:
        logging.error("Unsupported platform")
        return False

    if not destination.exists():
        if not download_file(url, destination):
            return False

    if extract_tar_gz(destination, directory):
        install_dir = directory / "MGLTools"
        install_script = directory / dir_name / "install.sh"
        cmd_args = ["sh", str(install_script), "-d", str(install_dir)]

        try:
            subprocess.check_output(
                cmd_args,
                cwd=directory / dir_name,
                shell=False,
                text=True,
                stderr=subprocess.DEVNULL,
            )
            logging.info("MGLTools installed successfully.")
            return True
        except subprocess.CalledProcessError as err:
            logging.error(f"Installation failed: {err}")
            return False
    else:
        return False
