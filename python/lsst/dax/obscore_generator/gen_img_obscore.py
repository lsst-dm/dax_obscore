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

import json
import math
import os
import urllib.parse as UrlParser
import uuid

import click
import lsst.geom as geom
import pandas as pd
from lsst.daf.base.dateTime.dateTime import DateTime
from lsst.daf.butler import Butler
from lsst.obs.base import Instrument

ROOT = os.path.dirname(__file__)
config_file = os.path.join(ROOT + "/config/gen_obscore_config.json")
with open(config_file) as f:
    _config = json.load(f)

_template = {
    "ds_bg3_repo": _config["ds_bg3_repo"],
    "dataproduct_type": _config["dataproduct_type"],
    "obs_collection": _config["obs_collection"],
    "facility_name": _config["facility_name"],
    "instrument_name": _config["instrument_name"],
}

_butler = Butler(_config["ds_bg3_repo"], collections=_config["obs_collection"])
_instrument = Instrument.fromName(_config["instrument_name"], _butler.registry)


"""
pr = cProfile.Profile()
pr.enable()
"""


@click.command()
@click.option("--ds_types")
@click.option("--out_dir")
def exec_gen_obscore(ds_types, out_dir):
    """Generate image obscore data in CSV."""
    print(f"Begin ObsCore data extraction from collection={_config['obs_collection']} ...")
    ds_type_list = ds_types.split(",")
    for dsType in ds_type_list:
        ds_collection = _config["obs_collection"]
        ds_obj = _butler.registry.queryDatasets(dsType, collections=[ds_collection])
        n = len(list(ds_obj))
        print(f"Processing ds={_config['ds_bg3_repo']} dsType={dsType} n={n} ...")
        df = pd.DataFrame([])
        i = 0
        for ref in _butler.registry.queryDatasets(dsType, collections=[ds_collection]):
            i += 1
            print(f"\rProcessing {i} of {n}", end="")
            r = _template.copy()
            r["core_id"] = str(uuid.uuid4())  # PRIMARY KEY
            r["planeID"] = str(uuid.uuid4())
            r["dataproduct_subtype"] = "lsst." + dsType
            r["calib_level"] = _config["calib_level_map"][dsType]
            r["target_name"] = None
            if dsType == "calexp":
                r["obs_id"] = ref.dataId["visit"]
                visit = _butler.registry.queryDimensionRecords("visit", dataId=ref.dataId)
                for d in visit:
                    r["target_name"] = d.target_name  # scheduler field ID
            elif dsType == "raw":
                r["obs_id"] = ref.dataId["exposure"]
                visit = _butler.registry.queryDimensionRecords("visit", dataId=ref.dataId)
                for d in visit:
                    r["target_name"] = d.target_name  # scheduler field ID
            elif dsType == "deepCoadd":
                r["obs_id"] = r["obs_collection"] + "-" + str(ref.dataId["patch"])
                r["target_name"] = str(ref.dataId["patch"])
            dataid = UrlParser.urlencode(ref.dataId)
            r["obs_publisher_did"] = "ivo://" + ds_collection + "/" + dsType + "?" + dataid
            uri = _butler.datastore.getURI(ref)
            r["access_url"] = uri
            exp = None
            try:
                wcs = _butler.get(dsType + ".wcs", ref.dataId)
            except (KeyError, LookupError):
                # note: costly to load the exposure below
                exp = _butler.getDirect(ref)
                wcs = exp.getWcs()
            try:
                bbox = _butler.get(dsType + ".bbox", ref.dataId)
                if bbox:
                    imageBox = geom.Box2D(bbox)
                else:
                    raise LookupError
            except (KeyError, LookupError):
                if not exp:
                    # note: costly to load the exposure below
                    exp = _butler.getDirect(ref)
                imageBox = geom.Box2D(exp.getBBox())
            visitInfo = None
            try:
                visitInfo = _butler.get(dsType + ".visitInfo", ref.dataId)
            except (KeyError, LookupError):
                if not exp:
                    # note: costly to load the exposure below
                    exp = _butler.getDirect(ref)
                expInfo = exp.getInfo()
                visitInfo = expInfo.getVisitInfo()
            if ref.dataId.hasRecords():
                r["lastmodified"] = str(ref.dataId.timespan.begin)
            else:
                r["lastmodified"] = DateTime.now().toString(DateTime.Timescale.UTC)
            try:
                f = _butler.get(dsType + ".filter", ref.dataId)
                fp = f.getFilterProperty()
                l_eff = fp.getLambdaEff()
                em_min = fp.getLambdaMin()
                if math.isnan(em_min):
                    em_min = l_eff
                em_max = fp.getLambdaMax()
                if math.isnan(em_max):
                    em_max = l_eff
                r["em_min"] = em_min
                r["em_max"] = em_max
                r["em_filter_name"] = ref.dataId["band"]
            except KeyError:
                print(f"band not found in dataId{ref.dataId}")
            # get uniform date and time at middle of exposure
            et = visitInfo.getExposureTime()
            if math.isnan(et) or et is None:
                r["t_exptime"] = et
            else:
                r["t_exptime"] = int(et)
            if ref.dataId.hasRecords():
                r["t_min"] = ref.dataId.timespan.begin.mjd
                r["t_max"] = ref.dataId.timespan.end.mjd
            else:
                r["t_min"] = None
                r["t_max"] = None
            corners = [wcs.pixelToSky(pt) for pt in imageBox.getCorners()]
            poly = ""
            for corner in corners:
                s_corner = str(corner)
                s_corner = s_corner.replace(",", "d,").replace(")", "d)")
                poly += s_corner
            poly = poly.replace(")(", "),(")
            r["position_bounds_spoly"] = "{" + poly + "}"
            poly = poly.strip("()")
            poly = poly.replace("),(", " ").replace(", ", " ")
            poly = poly.replace("d", "")
            r["s_region"] = "POLYGON ICRS " + poly
            r["s_pixel_scale"] = wcs.getPixelScale().asArcseconds()
            imageCenter = wcs.pixelToSky(imageBox.getCenter())
            r["s_ra"] = imageCenter.getRa().asRadians()
            r["s_dec"] = imageCenter.getDec().asRadians()
            w = imageBox.getWidth()
            h = imageBox.getHeight()
            radius = math.sqrt((w / 2) ** 2 + (h / 2) ** 2)
            # x_fov: diameter of bounding circle around s_ra, s_dec
            r["s_fov"] = 2 * radius * wcs.getPixelScale().asDegrees()
            # number of pixels in horizontal direction
            r["s_xel1"] = int(w)
            # number of pixels in vertical direction
            r["s_xel2"] = int(h)
            df = df.append(r, ignore_index=True)
        # fix the float output issue
        df["calib_level"] = df["calib_level"].astype(pd.Int64Dtype())
        df["t_exptime"] = df["t_exptime"].astype(pd.Int64Dtype())
        df["s_xel1"] = df["s_xel1"].astype(pd.Int64Dtype())
        df["s_xel2"] = df["s_xel2"].astype(pd.Int64Dtype())
        fn = os.path.join(out_dir, f"gen_obscore_{dsType}.csv")
        df.to_csv(fn, index=False, header=True)
        print("")  # finalizing newline
    print("ObsCore data generation complete")


"""
pr.disable()
pr.print_stats(sort="calls")
"""


if __name__ == "__main__":
    exec_gen_obscore()
