from osgeo import gdal, ogr, osr
from warnings import warn
import logging
import numpy as np
import os


logging.captureWarnings(True)
gdal.UseExceptions()


# TODO: Integrate some of the functionality from rasterio and/or shapely and/or fiona
class RavenGIS(object):
    """
    Base object for GIS analyses. Contains static methods that are used to perform data extractions of raster
    and vector data objects as well as performs comparison analyses between them.
    """

    def __init__(self):
        self._shape = None
        self._raster = None
        self._pretty_crs = None
        self._epsg = None

    def crs(self):
        return self._epsg

    @staticmethod
    def raster_bbox_corners(raster):
        """Returns the lat,lon boundaries for union of all features"""
        upper_left_x, x_res, x_skew, upper_left_y, y_skew, y_res = raster.GetGeoTransform()
        lower_right_x = upper_left_x + (raster.RasterXSize * x_res)
        lower_right_y = upper_left_y + (raster.RasterYSize * y_res)
        return upper_left_x, upper_left_y, lower_right_x, lower_right_y

    @staticmethod
    def shape_bbox_corners(bbox, transform):
        """Extracts the corner coordinates from a Geometry Reference"""
        xo = int(round((bbox[0] - transform[0]) / transform[1]))
        yo = int(round((transform[3] - bbox[3]) / transform[1]))
        xd = int(round((bbox[1] - bbox[0]) / transform[1]))
        yd = int(round((bbox[3] - bbox[2]) / transform[1]))
        return xo, yo, xd, yd

    @staticmethod
    def extract_crs(shape=None, raster=None):
        """Boilerplate for quickly parsing the CRS of given shape/raster layer."""
        if shape and raster is None:
            try:
                ogr_shape = ogr.Open(shape)
                assert ogr_shape is not None
            except AssertionError:
                raise Exception('Not a valid shape!')
            layer = ogr_shape.GetLayer()
            spa_ref = layer.GetSpatialRef()

            srs = osr.SpatialReference()
            srs.ImportFromWkt(spa_ref.ExportToWkt())

        elif raster and shape is None:
            try:
                gdal_raster = gdal.Open(raster, gdal.GA_ReadOnly)
                assert gdal_raster is not None
            except AssertionError:
                raise Exception('Not a valid raster!')
            srs = osr.SpatialReference()
            srs.ImportFromWkt(gdal_raster.GetProjection())

        else:
            raise Exception('Please provide a shape or a raster.')

        if srs.IsGeographic():
            structure = 'GEOGCS|AUTHORITY'
        else:
            structure = 'PROJCS|GEOGCS|AUTHORITY'

        epsg = srs.GetAttrValue(structure, 1)

        return srs.ExportToPrettyWkt(), epsg

    @staticmethod
    def open_raster(raster):
        """Tries to open a raster with GDAL"""
        if raster is None:
            raise Exception('No raster specified.')
        if isinstance(raster, str):
            if raster.endswith('.tif') or raster.endswith('.tiff'):
                gdal.GetDriverByName('GTiff')
        try:
            return gdal.Open(raster, gdal.GA_ReadOnly)
        except Exception as e:
            msg = "Failed to open DEM/raster image: {}".format(e)
            raise Exception(msg)

    @staticmethod
    def box(corners):
        """Creates a boundary box based on the extent of a raster"""
        ulx, uly, lrx, lry = corners
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(ulx, lry)
        ring.AddPoint(ulx, uly)
        ring.AddPoint(lrx, uly)
        ring.AddPoint(lrx, lry)
        ring.AddPoint(ulx, lry)
        box = ogr.Geometry(ogr.wkbPolygon)
        box.AddGeometry(ring)
        return box

    @staticmethod
    def open_shape(shape):
        if shape is None:
            raise Exception('No shape specified.')
        try:
            return ogr.Open(shape)
        except Exception as e:
            msg = 'Failed to load Shapefile/GML: {}'.format(e)
            raise Exception(msg)

    @staticmethod
    def _check_crs(shape=None, raster=None, crs=None):
        """
        Ensures that the CRS of the vector/raster and the specified CRS from EPSG string or shape/dem are compatible.
        Needed to circumvent problems that can arise from OGC-compliant WKT and ESRI-WKT CRS strings.
        Not sure of a better way to handle this issue.
        """

        try:
            if shape is None and raster is None:
                raise Exception("Please specify a shape or raster instance.")

            elif all([shape, raster, crs]):
                if len({shape, raster, crs}) == 1:
                    return True
                else:
                    return False

            elif (shape and raster and crs is None) or (shape and crs and raster is None)\
                    or (raster and crs and shape is None):

                if shape == raster:
                    return True
                elif shape == crs:
                    return True
                elif raster == crs:
                    return True
                else:
                    return False

            else:
                raise Exception("Please provide more arguments.")

        except Exception as e:
            msg = 'Failed to read SHAPE/DEM/CRS: {}'.format(e)
            raise Exception(msg)


class RavenRaster(RavenGIS):
    """
    Subclass of RavenGIS object

    A GIS object initialized with a raster/DEM. Prevents coordinate transformations of the DEM/raster
    which is important for discrete value categories (e.g. Land Use and Land Use Change data).
    """

    def __init__(self, raster):
        if isinstance(raster, str):
            super(RavenRaster, self).__init__()

            self._raster = self.open_raster(raster)
            self._filename = os.path.basename(raster)
            self._pretty_crs, self._epsg = self.extract_crs(raster=raster)
        else:
            raise Exception('File cannot be read!')

    def __repr__(self):
        return 'RavenGIS({!r}, {!r})'.format(self.__class__.__name__, self._filename)

    def __str__(self):
        return 'RavenGIS Raster: {}'.format(self._filename)

    def stats(self):
        """Outputs stats based on raster contents."""
        raster = self._raster
        transform = raster.GetGeoTransform()
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()

        msg = 'Raster stats: {}'.format(self._filename)
        print(msg)
        print('~' * len(msg))
        print('Layer boundary box: {}'.format(self.raster_bbox_corners(self._raster)))
        print('Geographic Transformation:{}'.format(transform), '\n',
              'Min value:{}'.format(array.min()), '\n',
              'Mean value:{}'.format(int(array.mean())), '\n',
              'Max value:{}'.format(array.max()))
        print('CRS:\n\t{}'.format(self._pretty_crs))

    def check_crs(self, shape=None, crs=None):
        """
        Ensures that the CRS of the vector/raster and the specified CRS from EPSG string or shape/dem are compatible.
        Needed to circumvent problems that can arise from OGC-compliant WKT and ESRI-WKT CRS strings.
        Not sure of a better way to handle this issue.
        """
        shape_crs = None
        if shape:
            shape_crs = self.extract_crs(shape=shape)[1]

        return self._check_crs(shape=shape_crs, raster=self._epsg, crs=crs)

    def _get_raster_counts(self):
        """A method for returning the cell values and frequencies of the RavenRaster"""
        raise NotImplementedError


class RavenShape(RavenGIS):
    """
    Subclass of RavenGIS object

    A GIS object initialized with a Shapefile/GML/GeoJSON. Prevents coordinate transformations of the shape
    which is important for calculating area and slope for specific regional projections.
    """

    # Examples via https://gis.stackexchange.com/a/195208/65343
    # Examples via https://stackoverflow.com/a/50039984/7322852
    # Examples via https://gis.stackexchange.com/a/126472/65343

    def __init__(self, shape):
        if isinstance(shape, str):
            super(RavenShape, self).__init__()

            self._shape = self.open_shape(shape)
            self._filename = os.path.basename(shape)
            self._pretty_crs, self._epsg = self.extract_crs(shape=shape)
        else:
            raise Exception('File cannot be read!')

    def __repr__(self):
        return 'RavenGIS({!r}, {!r})'.format(self.__class__.__name__, self._filename)

    def __str__(self):
        return 'RavenGIS Shape: {}'.format(self._filename)

    def __len__(self, *args, **kwargs):
        return self._get_feature_count()

    def stats(self):
        """Outputs stats based on shapefile contents."""
        msg = 'Shape stats: {}'.format(self._filename)
        print(msg)
        print('~' * len(msg))
        print('Layer name: {}'.format(self._get_layer_name()))
        print('Layer bounds: {}'.format(self._get_layer_bounds()))
        print('Number of features: {}'.format(self._get_feature_count()))
        print('Field names: {}'.format(self._get_fields()))
        print('CRS:\n\t{}'.format(self._pretty_crs))

    def _get_fields(self):
        """Returns the field names from attributes table of shape layer."""
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        fields = []
        for i in range(layer_def.GetFieldCount()):
            fields.append(layer_def.GetFieldDefn(i).GetName())
        return fields

    def _get_layer_name(self):
        """Returns the name of the shape layer. Can lead to surprising results with ESRI Shapefiles."""
        layer = self._shape.GetLayer()
        layer_def = layer.GetLayerDefn()
        return layer_def.GetName()

    def _get_layer_bounds(self):
        """Returns the lat,lon boundaries for union of all features in layer"""
        layer = self._shape.GetLayer()
        return layer.GetExtent()

    def _get_feature_count(self):
        """Returns count of all features in layer"""
        layer = self._shape.GetLayer()
        return layer.GetFeatureCount()

    def _get_feature_bounds(self):
        raise NotImplementedError

    def _get_layer_area(self):
        raise NotImplementedError

    def _get_feature_area(self):
        raise NotImplementedError

    def centroids(self):
        """Returns the centroids respecting all features of the shape layer."""
        layer = self._shape.GetLayer()

        for i, feature in enumerate(layer):
            geom = feature.GetGeometryRef()
            yield i, geom.Centroid().ExportToWkt()

        layer.ResetReading()

        return

    def check_crs(self, raster=None, crs=None):
        """
        Ensures that the CRS of the vector/raster and the specified CRS from EPSG string or shape/dem are compatible.
        Needed to circumvent problems that can arise from OGC-compliant WKT and ESRI-WKT CRS strings.
        Not sure of a better way to handle this issue.
        """
        raster_crs = None
        if raster:
            raster_crs = self.extract_crs(raster=raster)[1]

        return self._check_crs(shape=self._epsg, raster=raster_crs, crs=crs)

    def clip(self, raster):

        raster = self.open_raster(raster)
        transform = raster.GetGeoTransform()
        band = raster.GetRasterBand(1)

        # Create a shape to verify the spatial relationship between raster and shape
        raster_corners = self.raster_bbox_corners(raster)
        raster_bounds = self.box(raster_corners)

        shp = self._shape
        layer = shp.GetLayer()

        for i, feat in enumerate(layer):
            shape_geo = feat.GetGeometryRef()
            bbox = shape_geo.GetEnvelope()

            xo, yo, xd, yd = self.shape_bbox_corners(bbox, transform)
            arr = np.array([band.ReadAsArray(xo, yo, xd, yd)])

            if not raster_bounds.Intersects(shape_geo):
                msg = 'Clipped regions for feat {} do not intersect raster.'.format(i)
                warn(msg)
                continue
            elif not raster_bounds.Contains(shape_geo):
                msg = 'Clipped regions for feat {} are not contained by raster.'.format(i)
                warn(msg)
                continue

            yield i, arr

        layer.ResetReading()
        return

    def average(self, raster):
        clipper = self.clip(raster)

        # Create a dictionary with FIDs for clipped parcels
        for i, parcel in clipper:
            if isinstance(parcel, np.ndarray):
                yield i, np.mean(parcel)

        return
