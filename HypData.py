from typing import List
import pandas as pd
from HypTypes import *


np.seterr(divide='ignore', invalid='ignore')

""" Paths to all data files """
NMT_PATH_MAL = r'/home/rl/Projects/Eng/Hyperloop/NMT/malopolskie.txt'
NMT_PATH_MAZ = r'/home/rl/Projects/Eng/Hyperloop/NMT/mazowieckie.txt'
NMT_PATH_SWI = r'/home/rl/Projects/Eng/Hyperloop/NMT/swietokrzyskie.txt'
NMT_PATH_POD = r'/home/rl/Projects/Eng/Hyperloop/NMT/podkarpackie.txt'
NMT_PATH_LUB = r'/home/rl/Projects/Eng/Hyperloop/NMT/lubelskie.txt'
NMT_PATH_SLA = r'/home/rl/Projects/Eng/Hyperloop/NMT/slaskie.txt'
NMT_PATH_LOD = r'/home/rl/Projects/Eng/Hyperloop/NMT/lodzkie.txt'

FILTERED_DATA_PATH = r'FilteredData\\'
INVALID_DATA_PATH = r'InvalidData\\'

NMT_PATHS = [NMT_PATH_MAL, NMT_PATH_MAZ, NMT_PATH_SWI, NMT_PATH_POD, NMT_PATH_LUB, NMT_PATH_SLA, NMT_PATH_LOD]

REQ_POINTS = {
    DIST100M: ['000', '100', '200', '300', '400', '500', '600', '700', '800', '900'],
    DIST200M: ['000', '200', '400', '600', '800'],
    DIST500M: ['000', '500'],
    DIST1KM: ['000']
}

ALLOWED_DISTS = REQ_POINTS.keys()

""" Coordinates for Krakow and Warsaw """
X_KRK = 244000
Y_KRK = 567000
Z_KRK = 208
X_WAW = 487000
Y_WAW = 638000
Z_WAW = 108

""" Indexes for samples content """
YBEG = 3
YEND = 6
XBEG = 13
XEND = 16


def get_axis_xyz_data(x_start: int, y_start: int, range_km: int, dist_m: int) -> (np.array, np.array, np.array):
    """
    From raw data filter out not useful samples and convert to z values topology mesh
    :param x_start: x axis beginning point, must be in full kilometers
    :param y_start: y axis beginning point must be in full kilometers
    :param range_km: range in kilometers from which samples are requested
    :param dist_m: distance between fallowing points
    :return: array of samples in requested range and resolution in x,y,z dimentions
    """

    """ Check if parameters are correct """
    if x_start % DIST1KM != 0 or y_start % DIST1KM != 0 or dist_m not in ALLOWED_DISTS:
        raise ValueError('Incorrect parameter value!')

    """ Prefilter input data """
    nmt_data = []
    for nmt_path in NMT_PATHS:
        with open(nmt_path, 'r') as nmt_file:
            nmt_data += [sample for sample in nmt_file.readlines() if (sample[YBEG:YEND] in REQ_POINTS[dist_m] and
                                                                      sample[XBEG:XEND] in REQ_POINTS[dist_m])]
            print('Loaded samples from: {}'.format(nmt_path))
    print('Samples to be procced: {}'.format(len(nmt_data)))

    n_samples = np.int32(range_km * DIST1KM / dist_m)
    z_data = np.full((n_samples, n_samples), np.inf)
    print('Created empty mesh {}x{}'.format(n_samples, n_samples))

    x_req = [x_start + dist_m * n for n in range(n_samples)]
    y_req = [y_start + dist_m * n for n in range(n_samples)]
    print('Created requested indexes')

    i = 0
    for sample in nmt_data:
        y_sample, x_sample, z_sample = sample.split(' ')
        x_sample = np.int32(float(x_sample))
        y_sample = np.int32(float(y_sample))
        z_sample = np.int32(float(z_sample))
        i += 1
        if x_sample in x_req and y_sample in y_req:
            z_data[x_req.index(x_sample)][y_req.index(y_sample)] = z_sample
            print(int(i / len(nmt_data) * 100), '%\t', '|', x_req.index(x_sample), y_req.index(y_sample), '|', x_sample, y_sample, z_sample)

    return x_req, y_req, z_data


def get_axis_z_value(z_data: np.array, x_coordinate: int, y_coordinate: int, x_req: List, y_req: List) -> int:
    """
    Gets hight value for requested coordinates
    Params:                                                                     type:
    :param z_data: Axis x mesh                                                  np.array
    :param x_coordinate: x                                                      int
    :param y_coordinate: y                                                      int
    :param x_req: Axis x requested range                                        List
    :param y_req: Axis x requested range                                        List
    :return: Hight value
    """
    return z_data[list(x_req).index(round(x_coordinate))][list(y_req).index(round(y_coordinate))]


def read_z_data() -> np.array:
    z_data = np.array(pd.read_csv(FILTERED_DATA_PATH + 'z.csv'))
    z_data = np.array([elem[1:] for elem in z_data])
    return z_data


def read_filtered_data() -> (np.array, np.array, np.array):
    """
    Read preprocessed data
    Params:                                                                     type:
    :return: numpy array for every dimension                                    (np.array, np.array, np.array)
    """
    x_data = np.array(pd.read_csv(FILTERED_DATA_PATH + 'x.csv'))
    y_data = np.array(pd.read_csv(FILTERED_DATA_PATH + 'y.csv'))
    z_data = np.array(pd.read_csv(FILTERED_DATA_PATH + 'z.csv'))

    x_data = np.array([elem[1] for elem in x_data])
    y_data = np.array([elem[1] for elem in y_data])
    z_data = np.array([elem[1:] for elem in z_data])

    return x_data, y_data, z_data


def read_invalid_data() -> (np.array, np.array):
    x_inv = np.array(pd.read_csv(INVALID_DATA_PATH + 'x_inv.csv'))
    y_inv = np.array(pd.read_csv(INVALID_DATA_PATH + 'y_inv.csv'))

    x_inv = np.array([elem[1] for elem in x_inv])
    y_inv = np.array([elem[1] for elem in y_inv])

    return x_inv, y_inv

