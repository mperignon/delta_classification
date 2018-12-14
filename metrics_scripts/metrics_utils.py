import numpy as np
from osgeo import ogr, gdal, osr
from shapely.geometry import shape, Point, Polygon, LineString, MultiLineString, MultiPoint, MultiPolygon
from shapely.ops import polygonize_full
import cPickle as pickle

import fiona

def array2raster(filename,
                 rasterOrigin,
                 pixelWidth,
                 pixelHeight,
                 array,
                 spatial_ref = 32645):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]
    
    output_raster = (gdal.GetDriverByName('GTiff')
                         .Create(filename,
                                 cols, rows, 1,
                                 gdal.GDT_Float32))
                                                         
    output_raster.GetRasterBand(1).WriteArray( array ) 
    output_raster.SetGeoTransform((originX,
                                   pixelWidth, 0, originY, 0, pixelHeight))
                                   
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(spatial_ref)
    
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.FlushCache()
    
    
    
    
    
def latlong_to_index(pt, xmin, xmax, ymin, ymax, array_size):
    
    newx = int(array_size[1] * (pt[0] - xmin) / (xmax - xmin))
    newy = int(array_size[0] - array_size[0] * (pt[1] - ymin) / (ymax - ymin))
    
    return (newx, newy)
    
    
    
    
    
def create_shapefile_from_shapely_multi(features, filename,
                                        fields = {}, field_type = {},
                                        buffer_width = 0, spatial_ref = 32645):


    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(filename)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(spatial_ref)

    layer = ds.CreateLayer('', srs, ogr.wkbPolygon)


    for f in fields.keys():
        fieldDefn = ogr.FieldDefn(f, field_type[f])
        layer.CreateField(fieldDefn)


    defn = layer.GetLayerDefn()

    for i in range(len(features)):

        poly = features[i].buffer(buffer_width)
        feat = ogr.Feature(defn)

        for f in fields.keys():
            feat.SetField(f, fields[f][i])

        geom = ogr.CreateGeometryFromWkb(poly.wkb)
        feat.SetGeometry(geom)

        layer.CreateFeature(feat)
        feat = geom = None


    ds = layer = feat = geom = None
    



def create_tiff_from_shapefile(InputVector,
                               OutputImage,
                               RefImage,
                               burnValues = [1]):
    '''
    Turn a shapefile into a raster
    of the same shape and resolution as a RefImage
    with the values burnValues
    '''

    gdalformat = 'GTiff'
    datatype = gdal.GDT_Byte
    
    
    if not isinstance(burnValues, collections.Iterable) or isinstance(burnValues, (str, bytes)):
        burnValues = [burnValues]

    # Get projection info from reference image
    Image = gdal.Open(RefImage, gdal.GA_ReadOnly)

    # Open Shapefile
    Shapefile = ogr.Open(InputVector)
    Shapefile_layer = Shapefile.GetLayer()


    # Rasterise
    Output = gdal.GetDriverByName(gdalformat).Create(OutputImage,
                     Image.RasterXSize,
                     Image.RasterYSize,
                     1,
                     datatype,
                     options=['COMPRESS=DEFLATE'])
                     
                     
    Output.SetProjection(Image.GetProjectionRef())
    Output.SetGeoTransform(Image.GetGeoTransform()) 

    # Write data to band 1
    Band = Output.GetRasterBand(1)
    Band.SetNoDataValue(0)
    gdal.RasterizeLayer(Output, [1],
                        Shapefile_layer,
                        burn_values=burnValues)

    # Close datasets
    Band, Output, Image, Shapefile = None
    
    
    
def read_tiff_as_array(filename, get_info = True, normalize = False):
    
    src_ds = gdal.Open(filename)

    try:
        srcband = src_ds.GetRasterBand(1)
    except RuntimeError, e:
        # for example, try GetRasterBand(10)
        print 'Band ( %i ) not found' % band_num
        print e
        sys.exit(1)

    raster = srcband.ReadAsArray()
    
    if normalize:
        raster = (raster - np.min(raster)) / (np.max(raster) - np.min(raster))
    
    if get_info:
    
        ulx, xres, xskew, uly, yskew, yres  = src_ds.GetGeoTransform()
        lrx = ulx + (src_ds.RasterXSize * xres)
        lry = uly + (src_ds.RasterYSize * yres)
        
        srcband = src_ds = None
        
        return raster, (lrx, ulx, uly, lry, xres, yres)
    
    else:
    
        srcband = src_ds = None
        
        return raster
    
    
    
def get_value_of_raster_within_polygons(raster,
                                        polygons,
                                        r_xmin, r_xmax, 
                                        r_ymin, r_ymax, 
                                        min_valid_label = None):
    
    poly_points_utm = [i.representative_point().coords[0] for i in polygons]
    
    poly_points_index = ([latlong_to_index(i,
                          r_xmin, r_xmax, r_ymin, r_ymax, raster.shape)
                          for i in poly_points_utm])
    
    poly_points_value = [raster[i[::-1]] for i in poly_points_index]
    
    step = 10

    while (min_valid_label is not None) & (step < 25):

        bad_pts = np.where(np.array(poly_points_value) < min_valid_label)[0]

        if len(bad_pts) > 0:
            for i in bad_pts:

                pts = RegularGridSampling(polygons[i], step = step)
                pts_utm = [p.coords[0] for p in pts]
                
                pts_index = ([latlong_to_index(p, r_xmin, r_xmax, r_ymin,
                              r_ymax, raster.shape) for p in pts_utm])
                
                pts_value = [raster[p[::-1]] for p in pts_index]
                good_pts = np.where(np.array(pts_value) >= min_valid_label)[0]

                if len(good_pts) > 0:
                
                    poly_points_value[i] = pts_value[good_pts[0]]
                    poly_points_index[i] = pts_index[good_pts[0]]
                    
            step += 5
            
        else:
            step = 100
            
    poly_points_value = np.array(poly_points_value)  
            
    bad_pts = np.where(poly_points_value < min_valid_label)[0]
    
    return poly_points_value, poly_points_index, bad_pts
    
    
    
    

def RegularGridSampling(polygon,
                        x_interval = None,
                        y_interval = None,
                        step = None):
    """
    Perform sampling by substituting the polygon with a regular grid of
    sample points within it. The distance between the sample points is
    given by x_interval and y_interval.
    """
    
    samples = []
    
    if step is None:
        step = 10
    
    if x_interval is None:
        x_interval = (polygon.bounds[2] - polygon.bounds[0]) / step
        
    if y_interval is None:   
        y_interval = (polygon.bounds[3] - polygon.bounds[1]) / step
    
    ll = polygon.bounds[:2]
    ur = polygon.bounds[2:]
    
    low_x = int(ll[0]) / x_interval * x_interval
    upp_x = int(ur[0]) / x_interval * x_interval + x_interval
    
    low_y = int(ll[1]) / y_interval * y_interval
    upp_y = int(ur[1]) / y_interval * y_interval + y_interval

    for x in np.arange(low_x, upp_x, x_interval):
        for y in np.arange(low_y, upp_y, y_interval):
        
            p = Point(x, y)
            
            if p.within(polygon):
                samples.append(p)
                
    return MultiPoint(samples)





def SaveRaster(array,
                filename,
                r_xmin, r_ymin, res,
                downsampling = 1,
                spatial_ref = 32645):

    rasterOrigin = (r_xmin, r_ymin)
    pixelWidth = res * downsampling
    pixelHeight = res * downsampling

    reversed_arr = array[::-1]
    array2raster(filename,
                rasterOrigin,
                pixelWidth,
                pixelHeight,
                reversed_arr,
                spatial_ref = spatial_ref)
    
    
    
    
    
def azimuth(point1, point2):
    '''azimuth between 2 shapely points (interval 0 - 360), from vertical'''
    
    angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
    
    az = np.degrees(angle) if angle > 0 else np.degrees(angle) + 360
    
    return az    
    
    
    
def Polygon_axes(polygon):
    '''
    Calculates major and minor axes, and major axis azimuth,
    for a Shapely polygon
    '''

    rect = polygon.minimum_rotated_rectangle.exterior

    len1 = LineString(rect.coords[0:2]).length
    len2 = LineString(rect.coords[1:3]).length

    if len1 > len2:
        orientation = azimuth(Point(rect.coords[0]), Point(rect.coords[1]))
    else:
        orientation = azimuth(Point(rect.coords[1]), Point(rect.coords[2]))
        
    if orientation > 180:
        orientation = orientation - 180

    axminor, axmajor = np.sort([len1, len2])

    return axminor, axmajor, orientation    
    
    

def load_shapefile(filename, parameters = []):


    
    c = fiona.open(filename)

    if c[0]['geometry']['type'] == 'Polygon':
        shp = MultiPolygon([shape(pol['geometry']) for pol in c])

    elif c[0]['geometry']['type'] == 'LineString':
        shp = MultiLineString([shape(pol['geometry']) for pol in c])

    elif c[0]['geometry']['type'] == 'Point':
        shp = MultiPoint([shape(pol['geometry']) for pol in c])

    else:
        shp = [shape(pol['geometry']) for pol in c]
        
        
    
    if parameters is 'all':
        parameters = c[0]['properties'].keys()
        
    if type(parameters) is not list:
    
        parameters = list(parameters)


    shp_params = {}

    for param in parameters:
        shp_params[param] = [line['properties'][param] for line in c]

    c = None

    return shp, shp_params
    
    
    

def outline_to_mask(line, x, y):
    """Create mask from outline contour

    Parameters
    ----------
    line: array-like (N, 2)
    x, y: 1-D grid coordinates (input for meshgrid)

    Returns
    -------
    mask : 2-D boolean array (True inside)

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> poly = Point(0,0).buffer(1)
    >>> x = np.linspace(-5,5,100)
    >>> y = np.linspace(-5,5,100)
    >>> mask = outline_to_mask(poly.boundary, x, y)
    
    FROM: https://gist.github.com/perrette/a78f99b76aed54b6babf3597e0b331f8
    """
    
    import matplotlib.path as mplp
    
    mpath = mplp.Path(line)
    X, Y = np.meshgrid(x, y)
    points = np.array((X.flatten(), Y.flatten())).T
    mask = mpath.contains_points(points).reshape(X.shape)
    
    return mask



def load_zone_shapefiles(islands, files):
    
    zone_cat = np.zeros((len(islands),), dtype = 'int')
    
    for n,f in enumerate(files):
    
        c = fiona.open(f)
        zone_ids = [int(poly['properties']['id']) for poly in c]
        c = None

        zone_cat[np.array(zone_ids)] = n

    return zone_cat




def find_bounding_channels(shp,
                            poly,
                            save = True,
                            load_saved = False,
                            file_root = ''):

    if load_saved:
        bounds = pickle.load(open(
                        file_root + '/island_boundary_channels.p', "rb"))
                                  
        interior_channels = pickle.load(open(
                            file_root + '/island_interior_channels.p', "rb"))
                            
        channel_bounds = pickle.load(open(
                        file_root + '/island_channel_bounds.p', "rb"))
                                  
        flat_bounds = pickle.load(open(
                            file_root + '/island_flat_bounds.p', "rb"))



    else:

        midpts = [l.interpolate(0.5, normalized=True).buffer(5) for l in shp]

        bounds = []
        interior_channels = []

        # check if line midpoints intersect island outlines
        # to identify which lines make up each island
    
        for polygon in poly:
    
            touch = [i for i,l in enumerate(midpts)
                     if polygon.exterior.intersects(l)]
            bounds.append(touch)
        
            touch = [i for i,l in enumerate(midpts)
                     if polygon.contains(l)]
            interior_channels.append(touch)
            
            
        flat_bounds = np.unique([item for sublist in bounds for item in sublist])
        channel_bounds = MultiLineString([shp[i] for i in flat_bounds])

            
            
        if save:
        
            pickle.dump(bounds, open(
                        file_root + '/island_boundary_channels.p', "wb"))
                        
            pickle.dump(interior_channels, open(
                        file_root + '/island_interior_channels.p', "wb"))
                        
            pickle.dump(flat_bounds, open(
                        file_root + '/island_flat_bounds.p', "wb"))  
                        
            pickle.dump(interior_channels, open(
                        file_root + '/island_channel_bounds.p', "wb"))      
                        
            
        
    
    return bounds, interior_channels, flat_bounds, channel_bounds 
    
    
    
    
def find_nested_islands(poly, merge_islands = True):

    contained_islands = []

    for i in range(len(poly)):

        # if island has an inner ring
        if not poly[i].boundary.equals(poly[i].exterior):

            # check all other islands
            for j in range(len(poly)):
                if i != j:

                    # check if one is inside the other
                    inside = poly[j].within(Polygon(poly[i].exterior))

                    if inside:
                        contained_islands.append([i,j])

    
    if merge_islands:
            
        bad_islands = set([i[1] for i in contained_islands])
        holey_islands = set([i[0] for i in contained_islands])

        new_islands = [Polygon(poly[i].exterior)
                        if i in holey_islands else poly[i]
                        for i in range(len(poly))]
                        
        new_islands = [new_islands[i]
                        for i in range(len(new_islands))
                        if i not in bad_islands]

        poly = MultiPolygon(new_islands)
        
    
    return poly, contained_islands
    
    
    
    
    
    
def polygonize_shapely_lines(shp, size_threshold = 0):
    
    result, _,_,_ = polygonize_full(shp)
    polys = MultiPolygon([i for i in result if i.area > size_threshold])
    
    return polys
    
    
