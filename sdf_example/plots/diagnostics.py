import numpy as np
import os
from pylab import hist, plot as pyplot, xlabel, ylabel, xlim, ylim, savefig, acorr, mlab
from pylab import figure, subplot, subplots_adjust, gca, scatter, axvline, yticks, xticks


try:
    from statsmodels.regression.linear_model import yule_walker
    has_sm = True
except ImportError:
    has_sm = False


def spec(x, order=2):

    beta, sigma = yule_walker(x, order)
    return sigma**2 / (1. - np.sum(beta))**2



def geweke(x, first=.1, last=.5, intervals=20, maxlag=20):
    """Return z-scores for convergence diagnostics.

    Compare the mean of the first % of series with the mean of the last % of
    series. x is divided into a number of segments for which this difference is
    computed. If the series is converged, this score should oscillate between
    -1 and 1.

    Parameters
    ----------
    x : array-like
      The trace of some stochastic parameter.

    first : float
      The fraction of series at the beginning of the trace.

    last : float
      The fraction of series at the end to be compared with the section
      at the beginning.

    intervals : int
      The number of segments.

    maxlag : int
      Maximum autocorrelation lag for estimation of spectral variance

    Returns
    -------
    scores : list [[]]
      Return a list of [i, score], where i is the starting index for each
      interval and score the Geweke score on the interval.

    Notes
    -----
    The Geweke score on some series x is computed by:

      .. math:: \frac{E[x_s] - E[x_e]}{\sqrt{V[x_s] + V[x_e]}}

    where :math:`E` stands for the mean, :math:`V` the variance,

    :math:`x_s` a section at the start of the series and

    :math:`x_e` a section at the end of the series.

    References
    ----------
    Geweke (1992)
    """

    if not has_sm:
        print("statsmodels not available. Geweke diagnostic cannot be calculated.")
        return

    if np.ndim(x) > 1:
        return [geweke(y, first, last, intervals) for y in np.transpose(x)]

    # Filter out invalid intervals
    if first + last >= 1:
        raise ValueError(
            "Invalid intervals for Geweke convergence analysis", (first, last))

    # Initialize list of z-scores
    zscores = [None] * intervals

    # Starting points for calculations
    starts = np.linspace(0, int(len(x)*(1.-last)), intervals).astype(int)

    # Loop over start indices
    for i, s in enumerate(starts):

        # Size of remaining array
        x_trunc = x[s:]
        n = len(x_trunc)

        # Calculate slices
        first_slice = x_trunc[:int(first * n)]
        last_slice = x_trunc[int(last * n):]

        z = (first_slice.mean() - last_slice.mean())
        z /= np.sqrt(spec(first_slice)/len(first_slice) +
                     spec(last_slice)/len(last_slice))
        zscores[i] = len(x) - n, z

    return zscores





def geweke_plot(data, name, format='png', suffix='-diagnostic', path='./', fontmap=None):
    '''
    Generate Geweke (1992) diagnostic plots.

    Arguments:
        data: list
            List (or list of lists for vector-valued variables) of Geweke diagnostics, output
            from the `pymc.diagnostics.geweke` function .
        name: string
            The name of the plot.
        format (optional): string
            Graphic output format (defaults to png).
        suffix (optional): string
            Filename suffix (defaults to "-diagnostic").
        path (optional): string
            Specifies location for saving plots (defaults to local directory).
        fontmap (optional): dict
            Font map for plot.

    '''

    if fontmap is None:
        fontmap = {1: 10, 2: 8, 3: 6, 4: 5, 5: 4}

    # Generate new scatter plot
    figure()
    x, y = np.transpose(data)
    scatter(x.tolist(), y.tolist())

    # Plot options
    xlabel('First iteration', fontsize='x-small')
    ylabel('Z-score for %s' % name, fontsize='x-small')

    # Plot lines at +/- 2 sd from zero
    pyplot((np.min(x), np.max(x)), (2, 2), '--')
    pyplot((np.min(x), np.max(x)), (-2, -2), '--')

    # Set plot bound
    ylim(min(-2.5, np.min(y)), max(2.5, np.max(y)))
    xlim(0, np.max(x))

    # Save to file
    if not os.path.exists(path):
        os.mkdir(path)
    if not path.endswith('/'):
        path += '/'
    savefig("%s%s%s.%s" % (path, name, suffix, format))
