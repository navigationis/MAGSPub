import numpy as np

def course(lat1, long1, lat2, long2):

    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    long1 = np.radians(long1)
    long2 = np.radians(long2)
    a = np.sin(long1-long2)*np.cos(lat2)
    b = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(long1-long2)
    tc = -np.arctan2(a, b)
    if tc < 0: tc += 2*np.pi
    return np.degrees(tc)
    
def destination(tc, d, lat1, long1, r=3440):
    cd   = np.cos(d/r)
    sd   = np.sqrt(1 - cd**2)
    tc = np.radians(tc)
    lat1 = np.radians(lat1)
    long1 = np.radians(long1)
    cost = np.cos(tc)
    sint = np.sqrt(1 - cost**2)
    slt1 = np.sin(lat1)
    clt1 = np.cos(lat1)
    slat = slt1*cd + clt1*sd*cost
    dlon = np.arctan2(sint*sd*clt1, cd - slt1*slat)
    lat2 = np.arcsin(slat)
    long2 = long1 + dlon
    return np.degrees(lat2), np.degrees(long2)
    
def sph2plat(lon: np.ndarray, lat: np.ndarray, N: int = 1145) -> tuple[np.ndarray, np.ndarray]:
    # 1145 x 1145 map to force 0.02 degree pixels at equator
    r = (90 - lat) * N / 90
    x = r * np.cos(np.deg2rad(lon))
    y = r * np.sin(np.deg2rad(lon))
    return x.astype(np.int32), y.astype(np.int32)

def plat2sph(ix: np.ndarray, iy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lat = np.sqrt(ix**2 + iy**2) * 90 / 1145
    lon = np.rad2deg(np.arctan2(iy, ix))
    return lon, lat

def distance(long1, lat1, long2, lat2, r=3444):
    """
    Returns the distance (in nm, by default) between two points on a sphere.

    Parameters
    ----------
    long1   : float or ndarray
        Longitude of start point
    lat1    : float or ndarray
        Latitude of start point
    long2   : float or ndarray
        Longitude of end point
    lat2    : float or ndarray
        Latitude of end point
    r       : float
        Radius of sphere (Earth radius default, in nm. Change to preferred distance units)

    """
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    long1 = np.deg2rad(long1)
    long2 = np.deg2rad(long2)
    
    return 2*r*np.arcsin(np.sqrt((np.sin((lat1-lat2)/2))**2 + 
                         np.cos(lat1)*np.cos(lat2)*(np.sin((long1-long2)/2))**2))

class TabulatedFieldModel:
    def __init__(self, filename : str = None):
        if filename is None:
            from os.path import join, dirname
            filename = join(dirname(__file__), 'data', 'field-3m.npy')
        self.m = np.load(filename)
        self.RG = 1145
        self.RE = 3440

    def grad(self, lon: float, lat: float) -> tuple[float, float]:
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
        page = 0
        if lat < 0: 
            lat = -lat
            page = 1
        r = (np.pi/2 - lat) * self.RG / (np.pi/2)
        ix = (r * np.cos(lon)).astype(np.int32)
        iy = (r * np.sin(lon)).astype(np.int32)
        gx, gy = np.gradient(self.m[page,ix-1:ix+2,iy-1:iy+2])
        gx = gx[1,1]
        gy = gy[1,1]
        glat = -2/np.pi*(self.RG*gx*np.cos(lon) + self.RG*gy*np.sin(lon))
        glon =  r*(-gx*np.sin(lon) + gy*np.cos(lon))
        return glon / self.RE, glat / self.RE
    
class MAGSim:
    def __init__(self, field_model, baseline=10, rng=None):
        self.baseline = baseline
        self.field_model = field_model
        if rng is None: rng = np.random.default_rng()
        self.rng = rng

    def estimator(self, lon: float, lat: float, heading: float, speed: float):
        """
        Velocity estimator.
        """
        field_gradient = np.array(self.field_model.grad(lon, lat))
        hdg_r = np.deg2rad(90 - heading)
        u = np.array((np.cos(hdg_r), np.sin(hdg_r)))
        s = np.abs(np.dot(u, field_gradient))
        self._s = s
        dt = 0.00095 + 0.0061 / s + 0.0014 / s**2
        t = self.baseline / speed
        dv = speed * dt / t 
        return self.rng.normal(loc=speed, scale=dv)

