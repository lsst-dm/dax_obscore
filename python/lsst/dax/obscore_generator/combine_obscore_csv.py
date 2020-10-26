#!/usr/bin/env python
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#
# This file is part of dax_imgserv.
#
# Developed for the LSST Data Management System.
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import click
import glob
import pandas as pd


@click.command()
@click.option("--out_dir")
def exec_combine_obscore_csv(out_dir):
    file_pattern = f"{out_dir}/gen_obscore_"
    src_filenames = [i for i in glob.glob(f"{file_pattern}*")]
    print(f"Processing source files: {src_filenames}")
    combined_csv_data = pd.concat([pd.read_csv(f, delimiter=',', encoding='UTF-8') for f in src_filenames])
    out_file = f"{out_dir}/gen_obscore_out.csv"
    combined_csv_data.to_csv(out_file, index=False, header=True)
    print(f"Output: {out_file} ")


if __name__ == '__main__':
    exec_combine_obscore_csv()
