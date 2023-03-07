#include "geodesy.h"

#ifdef __cplusplus
extern "C" {
#endif

const ellipsoid_t sk42 =
{
    6378245.0,
    6356863.0,
    0.00335233
};

const ellipsoid_t pz90 =
{
	6378136.0,
	6356751.3618,
	0.003352803
};

const ellipsoid_t wgs84 =
{
    6378137.0,
    6356752.31424518,
    1.0 / 298.2572235630
};

const ellipsoid_t pz9011 =
{
    6378136.0,
    6356751.3618,
    0.003352803
};

const ecefconv_t wgs84_to_pz9011 =
{
    -0.013, 0.106, 0.022,
    -1.115071e-8, 1.71624e-8, -2.041066e-8,
    1.0 - 0.008e-6
};

const ecefconv_t pz9011_to_wgs84 =
{
    0.013, -0.106, -0.022,
    1.115071e-8, -1.71624e-8, 2.041066e-8,
    1.0 + 0.008e-6
};

const ecefconv_t pz9011_to_pz90 =
{
    1.443, -0.156, -0.222,
    1.115071e-8, -1.71624e-8, 6.506684e-7,
    1.0 + 0.228e-6
};

const ecefconv_t pz90_to_pz9011 =
{
    -1.443, 0.156, 0.222,
    -1.115071e-8, 1.71624e-8, -6.506684e-7,
    1 - 0.228e-6
};

const ecefconv_t sk42_to_pz9011 =
{
    23.557, -140.844, -79.778,
    -1.115071e-8, -1.679685e-6, -3.850439e-6,
    1.0 - 0.228e-6
};

const ecefconv_t pz9011_to_sk42 =
{
    -23.557, 140.844, 79.778,
    1.115071e-8, 1.679685e-6, 3.850439e-6,
    1.0 + 0.228e-6
};

static __inline double sqr(double x)
{
    return x * x;
}

double get_radius_normal(double lat_radians, const ellipsoid_t *ell)
{
    double aa, bb;
    /* Compute normal radius of planetary body */
    aa = ell->a * ell->a;
    bb = ell->b * ell->b;
    return aa / sqrt(aa * sqr(cos(lat_radians)) + bb * sqr(sin(lat_radians)));
}

void geodetic2ecef(double lat, double lon, double alt, double *x, double *y, double *z, const ellipsoid_t *ell)
{
    double N;
    /* radius of curvature of the prime vertical section */
    N = get_radius_normal(lat, ell);
    /*Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates. */
    *x = (N + alt) * cos(lat) * cos(lon);
    *y = (N + alt) * cos(lat) * sin(lon);
    *z = (N * sqr(ell->b / ell->a) + alt) * sin(lat);
}

void ecef2geodetic(double x, double y, double z, double *lat, double *lon, double *alt, const ellipsoid_t *ell)
{
    double ea, eb, rad, rho, vnew, v, c;
    int i;
    /*
    convert ECEF (meters) to geodetic coordinates
    input
    -----
    x,y,z  [meters] target ECEF location [0,Infinity)
    ell    reference ellipsoid
    output
    ------
    lat,lon   (radians)
    alt  (meters)
    Algorithm is based on
    http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    This algorithm provides a converging solution to the latitude equation
    in terms of the parametric or reduced latitude form (v)
    This algorithm provides a uniform solution over all latitudes as it does
    not involve division by cos(phi) or sin(phi)
    */
    ea = ell->a;
    eb = ell->b;
    rad = hypot(x, y);
    /* Constant required for Latitude equation */
    rho = atan2(eb * z, ea * rad);
    /* Constant required for latitude equation */
    c = (sqr(ea) - sqr(eb)) / hypot(ea * rad, eb * z);
    /* Starter for the Newtons Iteration Method */
    vnew = atan2(ea * z, eb * rad);
    /* Initializing the parametric latitude */
    v = 0;
    for (i = 0; i <= 5; i++)
    {
        v = vnew;
        /* Newtons Method for computing iterations */
        vnew = v - ((2 * sin(v - rho) - c * sin(2 * v)) /
                   (2 * (cos(v - rho) - c * cos(2 * v))));
        /* Computing latitude from the root of the latitude equation */
        *lat = atan2(ea * tan(vnew), eb);
        /* by inspection */
        *lon = atan2(y, x);
        *alt = (((rad - ea * cos(vnew)) * cos(*lat)) +
               ((z - eb * sin(vnew)) * sin(*lat)));
    }
}

void uvw2enu(double u, double v, double w, double lat0, double lon0, double *east, double *north, double *up)
{
    double t;
    t = cos(lon0) * u + sin(lon0) * v;
    *east = -sin(lon0) * u + cos(lon0) * v;
    *up = cos(lat0) * t + sin(lat0) * w;
    *north = -sin(lat0) * t + cos(lat0) * w;
}

void ecef2enu(double x, double y, double z, double lat0, double lon0, double h0, double *east, double *north, double *up, const ellipsoid_t *ell)
{
    double x0, y0, z0;
    /*
    input
    -----
    x,y,z  [meters] target ECEF location [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    output:
    -------
    e,n,u   East, North, Up [m]
    */
    geodetic2ecef(lat0, lon0, h0, &x0, &y0, &z0, ell);
    uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, east, north, up);
}

void geodetic2enu(double lat, double lon, double h, double lat0, double lon0, double h0, double *east, double *north, double *up, const ellipsoid_t *ell)
{
    double x1, y1, z1;
    double x2, y2, z2;
    /*
    input
    -----
    target: lat,lon, h
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    output:
    -------
    e,n,u   East, North, Up [m]
    */
    geodetic2ecef(lat, lon, h, &x1, &y1, &z1, ell);
    geodetic2ecef(lat0, lon0, h0, &x2, &y2, &z2, ell);
    uvw2enu(x1 - x2, y1 - y2, z1 - z2, lat0, lon0, east, north, up);
}

void enu2uvw(double east, double north, double up, double lat0, double lon0, double *u, double *v, double *w)
{
    double t;
    t = cos(lat0) * up - sin(lat0) * north;
    *w = sin(lat0) * up + cos(lat0) * north;
    *u = cos(lon0) * t - sin(lon0) * east;
    *v = sin(lon0) * t + cos(lon0) * east;
}

void enu2ecef(double e1, double n1, double u1, double lat0, double lon0, double h0, double *x, double *y, double *z, const ellipsoid_t *ell)
{
    double x0, y0, z0;
    double dx, dy, dz;
    /*
    Observer => Point
    inputs:
     e1, n1, u1 (meters)   east, north, up
     observer: lat0, lon0, h0 (degrees/radians,degrees/radians, meters)
    ell    reference ellipsoid
    output
    ------
    x,y,z  [meters] target ECEF location [0,Infinity)
    */
    geodetic2ecef(lat0, lon0, h0, &x0, &y0, &z0, ell);
    enu2uvw(e1, n1, u1, lat0, lon0, &dx, &dy, &dz);
    *x = x0 + dx;
    *y = y0 + dy;
    *z = z0 + dz;
}

void enu2geodetic(double e, double n, double up, double lat0, double lon0, double h0, double *lat, double *lon, double *h, const ellipsoid_t *ell)
{
    double x, y, z;
    /*
    input
    -----
    e,n,u   East, North, Up [m]
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    output:
    -------
    target: lat,lon, h  (degrees/radians,degrees/radians, meters)
    */
    enu2ecef(e, n, up, lat0, lon0, h0, &x, &y, &z, ell);
    ecef2geodetic(x, y, z, lat, lon, h, ell);
}

void geodetic2gk(double lat, double lon, int n, double *x, double *y)
{
    /*
    input
    -----
    target: lat,lon [radians]
    n: zone number
    output:
    -------
    projection: x,y (meters)
    */
    double sinB = sin(lat);
    double sin2B = sinB * sinB;
    double sin4B = sin2B * sin2B;
    double sin6B = sin4B * sin2B;
    double l = (lon * 180 / M_PI - (3 + 6 * (n - 1))) / 57.29577951;
    l = (l > M_PI) ? l - 2 * M_PI : l;
    l = (l <= -M_PI) ? l + 2 * M_PI : l;
    double l2 = l * l;
    *x = 6367558.4968 * lat - sin(2.0 * lat) * (16002.8900 + 66.9607 * sin2B + 0.3515 * sin4B
        - l2 * (1594561.25 +    5336.535 * sin2B + 26.790 * sin4B + 0.149 * sin6B
        + l2 * (672483.4 - 811219.9 * sin2B + 5420 * sin4B - 10.6 * sin6B
        + l2 * (278194 - 830174 * sin2B + 572434 * sin4B - 16010 * sin6B
        + l2 * (109500 - 574700 * sin2B + 863700 * sin4B - 398600 * sin6B)))));
    *y = 1e5 * (5 + 10 * n) + l * cos(lat) * (6378245 + 21346.1415 * sin2B + 107.1590 * sin4B + 0.5977 * sin6B
        + l2 * (1070204.16 - 2136826.66 * sin2B + 17.98 * sin4B - 11.99 * sin6B
        + l2 * (270806 - 1523417 * sin2B + 1327645 * sin4B - 21701 * sin6B
        + l2 * (79690 - 866190 * sin2B + 1730360 * sin4B - 945460 * sin6B))));
}

int geodetic2gk_zone(double lon)
{
    return (6.0 + lon * (180.0 / M_PI)) / 6.0;
}

void gk2geodetic(double x, double y, double *lat, double *lon)
{
	gkz2geodetic(x, y, gk_zone(y), lat, lon);
}

void gkz2geodetic(double x, double y, int n, double *lat, double *lon)
{
    double bet = x / 6367558.4968;
    double sin_bet = sin(bet);
    double sin2bet = sin_bet * sin_bet;
    double sin4bet = sin2bet * sin2bet;
    double B_0 = bet + sin(2.0 * bet) * (0.00252588685 - 0.00001491860 * sin2bet + 0.00000011904 * sin4bet);
    double z_0 = (y - 1e5 * (10 * n + 5)) / (6378245 * cos(B_0));
    double sin_B_0 = sin (B_0);
    double sin2B_0 = sin_B_0 * sin_B_0;
    double sin4B_0 = sin2B_0 * sin2B_0;
    double sin6B_0 = sin4B_0 * sin2B_0;
    double z_02 = z_0 * z_0;

    double delta_B = - z_02 * sin(2.0 * B_0) * (0.251684631 - 0.003369263 * sin2B_0 + 0.00001127 * sin4B_0
        - z_02 * (0.10500614 - 0.0455991 * sin2B_0 + 0.00228901 * sin4B_0 - 0.00002987 * sin6B_0
        - z_02 * (0.042858 - 0.025318 * sin2B_0 + 0.014346 *sin4B_0 - 0.001264 * sin6B_0
        - z_02 * (0.01672 - 0.00630 * sin2B_0 + 0.01188 * sin4B_0 - 0.00328 * sin6B_0))));
    *lat = B_0 + delta_B;

    double l = z_0 * (1.0 - 0.0033467108 * sin2B_0 - 0.0000056002 * sin4B_0 - 0.00000001897 * sin6B_0
        - z_02 * (0.16778975 + 0.16273586 * sin2B_0 - 0.00052490 * sin4B_0 - 0.00000846 * sin6B_0
        - z_02 * (0.0420025 + 0.1487407*sin2B_0 + 0.0059420 * sin4B_0 - 0.0000150 * sin6B_0
        - z_02 * (0.01225 + 0.09477 * sin2B_0 + 0.03282 * sin4B_0 - 0.00034 * sin6B_0
        - z_02 * (0.0038 + 0.0524 * sin2B_0 + 0.0482 * sin4B_0 - 0.0032 * sin6B_0)))));
    *lon = 6.0 * ((double)n - 0.5) / 57.29577951 + l;
}

int gk_zone(double y)
{
    // нет защиты от отрицательного y
    return (int)(y / 1e6);
}

void ecef2ecef(double src_x, double src_y, double src_z, double *dst_x, double *dst_y, double *dst_z, const ecefconv_t *conv)
{
    *dst_x = conv->m * (src_x + conv->wz * src_y - conv->wy * src_z ) + conv->dx;
    *dst_y = conv->m * (-src_x * conv->wz + src_y + conv->wx * src_z ) + conv->dy;
    *dst_z = conv->m * (src_x * conv->wy - conv->wx * src_y + src_z ) + conv->dz;
}

#ifdef __cplusplus
}
#endif
