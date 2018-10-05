import numpy as np
import matplotlib.pyplot as plt
from corner import corner

def flatten_chain(sampler, burnin=500, tailtrim=-1, acceptance_threshold=0.23):
    """
    Flatten sampler chains, excluding burnin period and chains below acceptrance threshold.
    """
    ndim = sampler.chain.shape[-1]
    return sampler.chain[sampler.acceptance_fraction >= acceptance_threshold, burnin:tailtrim, :].reshape((-1, ndim))

def plot_walkers(sampler, burnin=500, acceptance_threshold=0.23, labels=None):
    """
    Plot all walkers from sampler.
    """
    print('{:} chains below acceptance threshold ({:.2f})'.format(sum(sampler.acceptance_fraction <= acceptance_threshold),
                                                                  acceptance_threshold))

    ndim = sampler.chain.shape[-1]
    walkers = sampler.chain.shape[0]

    nplots = ndim + 1
    fig, axs = plt.subplots(nplots, figsize=(10, nplots * 1.25), sharex=True)
    if not isinstance(axs, np.ndarray):
        axs = [axs]

    if labels is None:
        labels = ['p{:.0f}'.format(i) for i in range(ndim)]
        
    for i in range(ndim):
        ax = axs[i]
        for walker in range(walkers):
            if sampler.acceptance_fraction[walker] >= acceptance_threshold:
                c='k'
                alpha_d = 4.
            else:
                c='r'
                alpha_d = 2.
            ax.plot(sampler.chain[walker][:, i], c=c, alpha=alpha_d / walkers)
        ax.set_ylabel(labels[i])
    
    ax = axs[-1]
    lnprob_settle = sampler.lnprobability[:, burnin:].ravel()
    lnprob_settle = np.nanmean(lnprob_settle[np.isfinite(lnprob_settle)])
    print(lnprob_settle)
    finite_prob = sampler.lnprobability.copy()
    finite_prob[~np.isfinite(finite_prob)] = np.nan
    xmin = burnin
    
    for walker in range(walkers):
        if sampler.acceptance_fraction[walker] >= acceptance_threshold:
            c='k'
            alpha_d = 4.
        else:
            c='r'
            alpha_d = 2.
        ax.plot(np.arange(xmin, sampler.lnprobability.shape[1]),
                sampler.lnprobability[walker, xmin:], c=c, alpha=alpha_d / walkers)

    ax.set_ylabel('lnprob')
    ax.set_ylim(1.2 * np.nanmin(finite_prob[sampler.acceptance_fraction >= acceptance_threshold, xmin:].ravel()), ax.get_ylim()[1])

    for ax in axs:
        ax.axvspan(0, burnin, alpha=0.1, color='b', zorder=-1)
        ax.set_xlim(0, sampler.chain.shape[1])
    
    axs[-1].set_xlabel('sampler step')
    fig.tight_layout()

    return fig, axs

def plot_corner(sampler, burnin=500, acceptance_threshold=0.23, bins=20, truths=None, labels=None, **kwargs):

    flatchain = flatten_chain(sampler, burnin=burnin, acceptance_threshold=acceptance_threshold)
    
    if labels is None:
        labels = ['p{:.0f}'.format(i) for i in range(len(flatchain))]

    fig = corner(flatchain, labels=labels, bins=bins, truths=truths, truth_color='r', **kwargs)
    
    return fig, fig.axes

def percentiles(sampler, percentiles=(2.5, 50, 97.5), burnin=500, acceptance_threshold=0.23):
    flatchain = flatten_chain(sampler, burnin=burnin, acceptance_threshold=acceptance_threshold)
    return np.percentile(flatchain, percentiles, 0).T