#!/usr/bin/env python-sirius

from fieldmaptrack import hallprobe as hall

magnets = ['BC-02',
           'BC-03',
           'BC-04',
           'BC-05',
           'BC-06',
           'BC-07',
           'BC-08',
           'BC-09',
           'BC-10',
           'BC-11',
           'BC-12',
           'BC-13',
           'BC-14',
           'BC-15',
           'BC-16',
           'BC-17',
           'BC-18',
           'BC-20',
           'BC-21']

fmap = ['2019-05-21_BC-02_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1744.dat',
        '2019-04-05_BC-03_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1314.dat',
        '2019-04-24_BC-04_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1472.dat',
        '2019-04-10_BC-05_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1364.dat',
        '2019-04-25_BC-06_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1489.dat',
        '2019-04-11_BC-07_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1379.dat',
        '2019-04-26_BC-08_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1506.dat',
        '2019-04-12_BC-09_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1396.dat',
        '2019-04-30_BC-10_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1556.dat',
        '2019-04-20_BC-11_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1446.dat',
        '2019-05-02_BC-12_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1575.dat',
        '2019-05-05_BC-13_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1602.dat',
        '2019-05-07_BC-14_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1617.dat',
        '2019-05-08_BC-15_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1631.dat',
        '2019-05-10_BC-16_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1655.dat',
        '2019-05-13_BC-17_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1671.dat',
        '2019-05-14_BC-18_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1681.dat',
        '2019-05-17_BC-20_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1713.dat',
        '2019-05-20_BC-21_Hall_X=-72_12mm_Z=-1500_1500mm_ID=1731.dat']

bc_path_fmt = '/home/imas/repos/si-dipoles-bc/model-13/analysis/hallprobe/production/x0-0p079mm-reftraj/{}/M1/{}'

tn = 4.2966
# t1 = 4.29564736842105
# t1_off_neg = 4.29843157894737
# t2_off_neg = 4.29652105263158

t1_off_pos = 4.29283157894737
t2_off_pos = 4.29675263157895


for i in range(0, len(magnets)):
    f = hall.DoubleFMapAnalysis(magnet=magnets[i], fmap_fname=fmap[i])
    f.energy = 3.0 * (t1_off_pos / tn) * (t2_off_pos / tn)
    f.traj_init_rx = 79.0e-3  # POSITIVO!!!

    if f.traj_rk_s_step < 0:
        f.traj_rk_s_step *= -1.0

    bc_path = bc_path_fmt.format(magnets[i], 'z-positive')
    f.create_input_rawfield(bc_path, force=True)
    f.create_input_trajectory(bc_path, force=True)
    f.create_input_multipoles(bc_path, force=True)
    f.create_input_model(bc_path, force=True)

    if f.traj_rk_s_step > 0:
        f.traj_rk_s_step *= -1.0

    bc_path = bc_path_fmt.format(magnets[i], 'z-negative')
    f.create_input_rawfield(bc_path, force=True)
    f.create_input_trajectory(bc_path, force=True)
    f.create_input_multipoles(bc_path, force=True)
    f.create_input_model(bc_path, force=True)
