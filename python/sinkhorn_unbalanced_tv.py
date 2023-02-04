from __future__ import division
import warnings
import numpy as np
from scipy.special import logsumexp

# from .utils import unif, dist

def sinkhorn_knopp_unbalanced_tv(a, b, M, reg, reg_m, numItermax=1000,
                              stopThr=1e-6, verbose=False, log=False, if_tv = False,**kwargs):
    r"""
    Solve the entropic regularization unbalanced optimal transport problem and return the loss

    The function solves the following optimization problem:

    .. math::
        W = \min_\gamma <\gamma,M>_F + reg\cdot\Omega(\gamma) + \reg_m KL(\gamma 1, a) + \reg_m KL(\gamma^T 1, b)

        s.t.
             \gamma\geq 0
    where :

    - M is the (dim_a, dim_b) metric cost matrix
    - :math:`\Omega` is the entropic regularization term :math:`\Omega(\gamma)=\sum_{i,j} \gamma_{i,j}\log(\gamma_{i,j})`
    - a and b are source and target unbalanced distributions
    - KL is the Kullback-Leibler divergence

    The algorithm used for solving the problem is the generalized Sinkhorn-Knopp matrix scaling algorithm as proposed in [10, 23]_


    Parameters
    ----------
    a : np.ndarray (dim_a,)
        Unnormalized histogram of dimension dim_a
    b : np.ndarray (dim_b,) or np.ndarray (dim_b, n_hists)
        One or multiple unnormalized histograms of dimension dim_b
        If many, compute all the OT distances (a, b_i)
    M : np.ndarray (dim_a, dim_b)
        loss matrix
    reg : float
        Entropy regularization term > 0
    reg_m: float
        Marginal relaxation term > 0
    numItermax : int, optional
        Max number of iterations
    stopThr : float, optional
        Stop threshol on error (> 0)
    verbose : bool, optional
        Print information along iterations
    log : bool, optional
        record log if True


    Returns
    -------
    if n_hists == 1:
        gamma : (dim_a x dim_b) ndarray
            Optimal transportation matrix for the given parameters
        log : dict
            log dictionary returned only if `log` is `True`
    else:
        ot_distance : (n_hists,) ndarray
            the OT distance between `a` and each of the histograms `b_i`
        log : dict
            log dictionary returned only if `log` is `True`
    Examples
    --------

    >>> import ot
    >>> a=[.5, .5]
    >>> b=[.5, .5]
    >>> M=[[0., 1.],[1., 0.]]
    >>> ot.unbalanced.sinkhorn_knopp_unbalanced(a, b, M, 1., 1.)
    array([[0.51122823, 0.18807035],
           [0.18807035, 0.51122823]])

    References
    ----------

    .. [10] Chizat, L., Peyré, G., Schmitzer, B., & Vialard, F. X. (2016).
        Scaling algorithms for unbalanced transport problems. arXiv preprint
        arXiv:1607.05816.

    .. [25] Frogner C., Zhang C., Mobahi H., Araya-Polo M., Poggio T. :
        Learning with a Wasserstein Loss,  Advances in Neural Information
        Processing Systems (NIPS) 2015

    See Also
    --------
    ot.lp.emd : Unregularized OT
    ot.optim.cg : General regularized OT

    """

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)

    dim_a, dim_b = M.shape

    if len(a) == 0:
        a = np.ones(dim_a, dtype=np.float64) / dim_a
    if len(b) == 0:
        b = np.ones(dim_b, dtype=np.float64) / dim_b

    if len(b.shape) > 1:
        n_hists = b.shape[1]
    else:
        n_hists = 0

    if log:
        log = {'err': []}

    # we assume that no distances are null except those of the diagonal of
    # distances
    if n_hists:
        u = np.ones((dim_a, 1)) / dim_a
        v = np.ones((dim_b, n_hists)) / dim_b
        a = a.reshape(dim_a, 1)
    else:
        u = np.ones(dim_a) / dim_a
        v = np.ones(dim_b) / dim_b

    # Next 3 lines equivalent to K= np.exp(-M/reg), but faster to compute
    K = np.empty(M.shape, dtype=M.dtype)
    np.divide(M, -reg, out=K)
    np.exp(K, out=K)

    fi = reg_m / (reg_m + reg)

    err = 1.

    for i in range(numItermax):
        uprev = u
        vprev = v

        # change this part for the TV-div
        Kv = K.dot(v)
        if if_tv:
            u_potential = (a / Kv)
            u = clamp(u_potential, reg_m)
        else:
            u = (a / Kv) ** fi
        Ktu = K.T.dot(u)
        if if_tv:
            v_potential = (b / Ktu)
            v = clamp(v_potential, reg_m)
        else:
            v = (b / Ktu) ** fi

        if (np.any(Ktu == 0.)
                or np.any(np.isnan(u)) or np.any(np.isnan(v))
                or np.any(np.isinf(u)) or np.any(np.isinf(v))):
            # we have reached the machine precision
            # come back to previous solution and quit loop
            warnings.warn('Numerical errors at iteration %s' % i)
            u = uprev
            v = vprev
            break

        err_u = abs(u - uprev).max() / max(abs(u).max(), abs(uprev).max(), 1.)
        err_v = abs(v - vprev).max() / max(abs(v).max(), abs(vprev).max(), 1.)
        err = 0.5 * (err_u + err_v)
        if log:
            log['err'].append(err)
            if verbose:
                if i % 50 == 0:
                    print(
                        '{:5s}|{:12s}'.format('It.', 'Err') + '\n' + '-' * 19)
                print('{:5d}|{:8e}|'.format(i, err))
        if err < stopThr:
            break

    if log:
        log['logu'] = np.log(u + 1e-300)
        log['logv'] = np.log(v + 1e-300)

    if n_hists:  # return only loss
        res = np.einsum('ik,ij,jk,ij->k', u, K, v, M)
        if log:
            return res, log
        else:
            return res

    else:  # return OT matrix

        if log:
            return u[:, None] * K * v[None, :], log
        else:
            return u[:, None] * K * v[None, :]


def clamp(potential, reg_m):
    clamped = potential
    neg_idx = [i for i in range(len(potential)) if potential[i] < -reg_m]
    pos_idx = [i for i in range(len(potential)) if potential[i] > reg_m]
    clamped[neg_idx] = -reg_m
    clamped[pos_idx] = reg_m

    return clamped