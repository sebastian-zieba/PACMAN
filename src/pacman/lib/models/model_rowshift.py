import sys
sys.path.insert(0,'..')

def model_rowshift(t, data, params, visit):
    rowshift_vf, rowshift_vr = params
    rowshift_vf = rowshift_vf[visit]#slope for forward scans
    rowshift_vr = rowshift_vr[visit]#slope for reverse scans

    return 1 + rowshift_vr*data.rowshift[data.vis_idx[visit]]*data.scan_direction[data.vis_idx[visit]] \
             + rowshift_vf*data.rowshift[data.vis_idx[visit]]*(1-data.scan_direction[data.vis_idx[visit]])