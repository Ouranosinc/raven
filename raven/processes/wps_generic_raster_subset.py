import logging
import tempfile

from pywps import ComplexOutput
from pywps import LiteralInput, ComplexInput
from pywps import Process, FORMATS, Format
from pywps.app.Common import Metadata
from rasterstats import zonal_stats

from raven.utilities import gis
from raven.utils import (
    archive_sniffer,
    crs_sniffer,
    single_file_check,
    raster_datatype_sniffer,
    generic_raster_warp,
    zonalstats_raster_file,
)

LOGGER = logging.getLogger("PYWPS")


class RasterSubsetProcess(Process):
    """Given files containing vector data and raster data, perform zonal statistics of the overlapping regions."""

    def __init__(self):
        inputs = [
            ComplexInput(
                "shape",
                "Vector Shape",
                abstract="An ESRI Shapefile, GML, JSON, GeoJSON, or single layer GeoPackage."
                " The ESRI Shapefile must be zipped and contain the .shp, .shx, and .dbf."
                " The shape and raster should have a matching CRS.",
                min_occurs=1,
                max_occurs=1,
                supported_formats=[
                    FORMATS.GEOJSON,
                    FORMATS.GML,
                    FORMATS.JSON,
                    FORMATS.SHP,
                ],
            ),
            ComplexInput(
                "raster",
                "Gridded raster data set",
                abstract="The raster to be queried. Defaults to the EarthEnv-DEM90 product.",
                metadata=[
                    Metadata("EarthEnv-DEM90", "https://www.earthenv.org/DEM"),
                    Metadata(
                        "Robinson, Natalie, James Regetz, and Robert P. Guralnick (2014). "
                        "EarthEnv-DEM90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m Digital "
                        "Elevation Model from Fused ASTER and SRTM Data. ISPRS Journal of "
                        "Photogrammetry and Remote Sensing 87: 57–67.",
                        "https://doi.org/10.1016/j.isprsjprs.2013.11.002",
                    ),
                ],
                min_occurs=0,
                max_occurs=1,
                supported_formats=[FORMATS.GEOTIFF],
            ),
            LiteralInput(
                "band",
                "Raster band",
                data_type="integer",
                default=1,
                abstract="Band of raster examined to perform zonal statistics. Default: 1",
                min_occurs=1,
                max_occurs=1,
            ),
            LiteralInput(
                "select_all_touching",
                "Additionally select boundary pixels that are touched by shape",
                data_type="boolean",
                default="false",
            ),
        ]

        outputs = [
            ComplexOutput(
                "raster",
                "DEM subset of `shape` region in GeoTIFF format.",
                abstract="Elevation statistics: min, max, mean, median, sum, nodata",
                as_reference=True,
                supported_formats=[
                    Format("application/zip"),
                ],
            ),
        ]

        super(RasterSubsetProcess, self).__init__(
            self._handler,
            identifier="raster-subset",
            title="Raster Subset",
            version="1.0",
            abstract="Return a masked raster based on boundaries of a vector file.",
            metadata=[],
            inputs=inputs,
            outputs=outputs,
            status_supported=True,
            store_supported=True,
        )

    def _handler(self, request, response):

        shape_url = request.inputs["shape"][0].file
        band = request.inputs["band"][0].data
        touches = request.inputs["select_all_touching"][0].data

        vectors = [".gml", ".shp", ".gpkg", ".geojson", ".json"]
        vector_file = single_file_check(
            archive_sniffer(shape_url, working_dir=self.workdir, extensions=vectors)
        )

        if "raster" in request.inputs:
            raster_url = request.inputs["raster"][0].file
            rasters = [".tiff", ".tif"]
            raster_file = single_file_check(
                archive_sniffer(
                    raster_url, working_dir=self.workdir, extensions=rasters
                )
            )
        else:
            bbox = gis.get_bbox(vector_file)
            raster_url = "public:EarthEnv_DEM90_NorthAmerica"
            raster_bytes = gis.get_raster_wcs(bbox, geographic=True, layer=raster_url)
            raster_file = tempfile.NamedTemporaryFile(
                prefix="wcs_", suffix=".tiff", delete=False, dir=self.workdir
            ).name
            with open(raster_file, "wb") as f:
                f.write(raster_bytes)

        vec_crs, ras_crs = crs_sniffer(vector_file), crs_sniffer(raster_file)

        if ras_crs != vec_crs:
            msg = f"CRS for files {vector_file} and {raster_file} are not the same. Reprojecting raster..."
            LOGGER.warning(msg)

            projected = tempfile.NamedTemporaryFile(
                prefix="reprojected_", suffix=".json", delete=False, dir=self.workdir
            ).name
            generic_raster_warp(vector_file, projected, target_crs=vec_crs)
            raster_file = projected

        data_type = raster_datatype_sniffer(raster_file)

        try:
            stats = zonal_stats(
                vector_file,
                raster_file,
                band=band,
                all_touched=touches,
                raster_out=True,
            )

            response.outputs["raster"].file = zonalstats_raster_file(
                stats,
                working_dir=self.workdir,
                identifier=self.identifier,
                data_type=data_type,
                crs=vec_crs or ras_crs,
            )

            # for i in range(len(stats)):
            #
            #     file = 'subset_{}.tiff'.format(i + 1)
            #     raster_subset = os.path.join(out_dir, file)
            #
            #     try:
            #         raster_location = stats[i]
            #         raster = raster_location['mini_raster_array']
            #         grid_properties = raster_location['mini_raster_affine'][0:6]
            #         nodata = raster_location['mini_raster_nodata']
            #
            #         aff = Affine(*grid_properties)
            #
            #         LOGGER.info('Writing raster data to {}'.format(raster_subset))
            #
            #         masked_array = np.ma.masked_values(raster, nodata)
            #         if masked_array.mask.all():
            #             msg = 'Subset {} is empty, continuing...'.format(i)
            #             LOGGER.warning(msg)
            #
            #         normal_array = np.asarray(masked_array, dtype=data_type)
            #
            #         # Write to GeoTIFF
            #         with rio.open(raster_subset, 'w', driver='GTiff', count=1, compress=raster_compression,
            #                       height=raster.shape[0], width=raster.shape[1], dtype=data_type, transform=aff,
            #                       crs=vec_crs or ras_crs, nodata=nodata) as f:
            #             f.write(normal_array, 1)
            #
            #     except Exception as e:
            #         msg = 'Failed to write raster outputs: {}'.format(e)
            #         LOGGER.error(msg)
            #         raise Exception(msg)
            #
            # # `shutil.make_archive` could potentially cause problems with multi-thread? Worth investigating later.
            # out_fn = os.path.join(self.workdir, self.identifier)
            # shutil.make_archive(base_name=out_fn, format='zip', root_dir=out_dir, logger=LOGGER)
            #
            # response.outputs['raster'].file = '{}.zip'.format(out_fn)

        except Exception as e:
            msg = f"Failed to perform raster subset using {shape_url} and {raster_url}: {e}"
            LOGGER.error(msg)
            raise Exception(msg)

        return response
