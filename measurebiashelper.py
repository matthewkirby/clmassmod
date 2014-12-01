import stats

def LogSum2DGaussianWrapper(args):

    return stats.LogSum2DGaussian(xs = args['cluster_like_samples'],
                                      mu = args['cluster_mean'],
                                      invcovar = args['cluster_invcovar'],
                                      sqrtdetcovar = args['cluster_detcovar'])

###################
