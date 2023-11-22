import numpy as np


def constants_cj(t, data, params, visit = 0):
    """
    Example:
        In [47]: iexp_orb_sp = np.array([0,1,2,3,0,1,2,3,4])

        In [48]: Cs = np.array([[7.8], [8.3], [8.5], [8.6], [8.65]])

        In [49]: C_data_mask = [iexp_orb_sp == i for i in range(max(iexp_orb_sp)+1)]

        In [50]: C_data_mask
        Out[50]:
        [array([ True, False, False, False,  True, False, False, False, False]),
         array([False,  True, False, False, False,  True, False, False, False]),
         array([False, False,  True, False, False, False,  True, False, False]),
         array([False, False, False,  True, False, False, False,  True, False]),
         array([False, False, False, False, False, False, False, False,  True])]

        In [51]: C_data_mask*Cs
        Out[51]:
        array([[7.8 , 0.  , 0.  , 0.  , 7.8 , 0.  , 0.  , 0.  , 0.  ],
               [0.  , 8.3 , 0.  , 0.  , 0.  , 8.3 , 0.  , 0.  , 0.  ],
               [0.  , 0.  , 8.5 , 0.  , 0.  , 0.  , 8.5 , 0.  , 0.  ],
               [0.  , 0.  , 0.  , 8.6 , 0.  , 0.  , 0.  , 8.6 , 0.  ],
               [0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 8.65]])

        In [52]: np.sum(C_data_mask*Cs, axis=0)
        Out[52]: array([7.8 , 8.3 , 8.5 , 8.6 , 7.8 , 8.3 , 8.5 , 8.6 , 8.65])

    """
    cjs = np.zeros((len(params),1))
    for j in range(len(params)):
        cjs[j] = 10. ** (params[j][visit])

    C_data_mask = np.array([data.iexp_orb_sp == i for i in range(max(data.iexp_orb_sp)+1)])

    return np.sum(C_data_mask * cjs, axis=0)
