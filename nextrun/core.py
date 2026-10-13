"""Pure numerical and causal state functions; no formal seed generation here."""
from dataclasses import dataclass
from statistics import NormalDist
import numpy as np

Z = NormalDist().inv_cdf(0.975)
TRI = np.triu_indices(20)
A_METHODS = ('A', 'A_full', 'A_clip', 'A_order', 'A_inf')


def frozen(a, dtype=np.float64):
    a = np.array(a, dtype=dtype, copy=True, order='C')
    a.flags.writeable = False
    assert a.flags.owndata
    return a


def euler_step(state, time, h):
    x, u, v = state
    return np.array([x + h*(np.sin(time)-0.1*x+0.5*np.tanh(u)),
                     u + h*(x-0.2*u-2*np.pi*v),
                     v + h*(2*np.pi*u-0.2*v)], dtype=np.float64)


def trajectory(h, burn=500, count=3221):
    stride = round(0.1/h)
    assert abs(h*stride - 0.1) < 1e-15
    state = np.zeros(3)
    out = []
    for k in range((burn+count)*stride):
        state = euler_step(state, k*h, h)
        if (k+1) % stride == 0 and (k+1)//stride > burn:
            out.append(state.copy())
    out = np.array(out)
    assert out.shape == (count, 3) and np.isfinite(out).all()
    return out


def observe(states, innovation, sigma):
    return np.asarray(states)[..., 0] + np.asarray(sigma)*innovation


def indices(offset=21):
    target = np.arange(3200)+offset
    inputs = target[:, None]-np.arange(1,21)
    return inputs, target


def noise_scales(offset, condition):
    sigma = np.full(3221, 0.15)
    if condition == 'shift':
        sigma[2400+160+offset:2400+800+offset] = 0.35
    else:
        assert condition == 'constant'
    return sigma


def feature(lags):
    assert lags.shape[-1] == 20 and np.isfinite(lags).all()
    return np.concatenate([np.ones(lags.shape[:-1]+(1,)), lags,
                           lags[..., TRI[0]]*lags[..., TRI[1]]], axis=-1)


@dataclass(frozen=True)
class Model:
    w: np.ndarray
    chol: np.ndarray
    tau: float
    df: float
    tau_df: float


def fit(train_lags, targets):
    assert len(targets) == len(train_lags) == 1600
    r = feature(train_lags)
    chol = np.linalg.cholesky(r.T@r + 16*np.eye(231))
    w = np.linalg.solve(chol.T, np.linalg.solve(chol, r.T@targets))
    resid = targets-r@w
    tau = float(np.sqrt(np.mean(resid**2)))
    # Solve the prescribed triangular equations; no explicit matrix inverse.
    basis_solutions = np.linalg.solve(chol, np.eye(231))
    df = float(231-16*np.sum(basis_solutions**2))
    assert 0 <= df <= 231
    tau_df = float(np.sqrt(np.sum(resid**2)/(1600-df)))
    assert 1 <= tau_df/tau <= 40/37
    assert np.isfinite(w).all()
    return Model(frozen(w), frozen(chol), tau, df, tau_df)


def predict(model, past_lags):
    assert past_lags.flags.owndata and not past_lags.flags.writeable
    r = feature(past_lags)
    z = np.linalg.solve(model.chol, r)
    leverage = float(z@z)
    assert leverage >= 0
    return float(r@model.w), leverage


@dataclass(frozen=True)
class History:
    values: np.ndarray
    ids: np.ndarray
    alpha: float = 0.05


def history(values, ids, alpha=0.05, absolute=False):
    values = np.abs(values) if absolute else values
    if absolute:
        assert np.all(np.asarray(values) >= 0)
    return History(frozen(values), frozen(ids, np.int64), float(alpha))


def initial_history(calibration_resid, method):
    full = method == 'A_full'
    vals = calibration_resid if full else calibration_resid[-100:]
    ids = np.arange(-len(vals), 0)
    return history(vals, ids, absolute=method != 'G')


def quantile_width(vals, alpha, method):
    assert np.all(vals >= 0) and np.isfinite(vals).all()
    p = 1-alpha
    if not (0 <= alpha <= 1 and 0 <= p <= 1):
        raise ValueError('invalid quantile probability')
    if method == 'A_order':
        k = int(np.ceil((len(vals)+1)*p))
        if k == len(vals)+1:
            return np.inf
        if k == 0:
            return 0.0
        return float(np.partition(vals, k-1)[k-1])
    if method == 'A_inf' and alpha == 0:
        return np.inf
    return float(np.quantile(vals, p, method='linear'))


def halfwidth(method, model, leverage, fixed_q, state=None):
    if method == 'B':
        h = Z*model.tau*np.sqrt(1+leverage)
    elif method == 'B_df':
        h = Z*model.tau_df*np.sqrt(1+leverage)
    elif method == 'F':
        h = fixed_q
    elif method in ('G', 'G_L'):
        h = Z*np.sqrt(np.mean(state.values**2))
        if method == 'G_L':
            h *= np.sqrt(1+leverage)
    elif method == 'rolling_F':
        h = np.partition(state.values, 95)[95]
    else:
        h = quantile_width(state.values, state.alpha, method)
    if h < 0 or np.isnan(h) or (not np.isfinite(h) and method not in ('A_inf','A_order')):
        raise ValueError('invalid interval halfwidth')
    return float(h)


def reveal(state, method, residual, miss, t):
    assert t > state.ids[-1]
    absolute = method not in ('G','G_L')
    vals = np.append(state.values, abs(residual) if absolute else residual)
    ids = np.append(state.ids, t)
    if method != 'A_full':
        vals, ids = vals[-100:], ids[-100:]
    alpha = state.alpha
    if method in A_METHODS:
        low, high = (0.001,0.999) if method == 'A_clip' else (0,1)
        alpha = float(np.clip(alpha+0.01*(0.05-int(miss)), low, high))
    return history(vals, ids, alpha, absolute)


def score(y, lower, upper):
    if np.isnan(lower) or np.isnan(upper) or lower > upper:
        raise ValueError('invalid endpoints')
    hit = bool(lower <= y <= upper)
    if np.isneginf(lower) and np.isposinf(upper):
        return hit, np.inf, np.inf
    if not np.isfinite([y,lower,upper]).all():
        raise ValueError('unexpected non-finite endpoint')
    width = upper-lower
    return hit, width, width+40*max(lower-y, 0)+40*max(y-upper, 0)


def evaluate(model, y, target_ids, cal_resid, methods):
    """Only the orchestrator holds y. Prediction functions receive owning past copies.

    Each log row's center/endpoints/alpha/history IDs is fixed before y[q] is read.
    At t, target_ids[:t] have been revealed; lag indices are all strictly < q.
    """
    m = len(methods)
    log = {k: np.empty((m,800)) for k in
           ('lower','upper','halfwidth','hit','width','winkler','alpha','history_max','history_first','history_last','infinite')}
    log.update(center=np.empty(800), leverage=np.empty(800), target=np.empty(800))
    states = {method: initial_history(cal_resid, method) for method in methods if method not in ('B','F','B_df')}
    fixed_q = float(np.partition(np.abs(cal_resid), 760)[760])
    for t, q in enumerate(target_ids):
        lag_ids = q-np.arange(1,21)
        assert np.all(lag_ids < q)
        past = frozen(y[lag_ids])
        center, lev = predict(model, past)
        log['center'][t], log['leverage'][t] = center, lev
        for k, method in enumerate(methods):
            state = states.get(method)
            if method in ('G','A') and t >= 100:
                assert np.array_equal(state.ids, np.arange(t-100,t))
            h = halfwidth(method, model, lev, fixed_q, state)
            # The entire prediction entry is committed before the reveal below.
            log['lower'][k,t], log['upper'][k,t] = center-h, center+h
            log['halfwidth'][k,t], log['infinite'][k,t] = h, np.isinf(h)
            log['alpha'][k,t] = state.alpha if method in A_METHODS else np.nan
            log['history_max'][k,t] = np.max(np.abs(state.values)) if state else np.nan
            log['history_first'][k,t] = state.ids[0] if state else np.nan
            log['history_last'][k,t] = state.ids[-1] if state else np.nan
        current_y = float(y[q])  # reveal after all prediction entries are fixed
        log['target'][t] = current_y
        for k, method in enumerate(methods):
            hit, width, wink = score(current_y, log['lower'][k,t], log['upper'][k,t])
            log['hit'][k,t], log['width'][k,t], log['winkler'][k,t] = hit, width, wink
            if method in states:
                states[method] = reveal(states[method], method, current_y-center, not hit, t)
        if 'G_L' in methods:
            g, gl = methods.index('G'), methods.index('G_L')
            assert log['lower'][gl,t] <= log['lower'][g,t] and log['upper'][gl,t] >= log['upper'][g,t]
            assert log['hit'][gl,t] >= log['hit'][g,t]
    return log
