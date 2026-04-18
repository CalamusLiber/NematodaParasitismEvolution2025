#!/apps/miniforge3/bin/python3.12
#SBATCH --job-name=SumMCMCTr
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# ============================================================================
# SumMCMCtree.py
# Purpose: Summarize and visualize MCMCtree posterior and prior results.
#   - read MCMCtree stdout and MCMC sample files
#   - construct ArviZ inference data for diagnostics
#   - compare prior and posterior node age distributions
#   - build dated phylogenetic trees with branch comments
#   - export summary tables, plots, and Nexus tree files
# ============================================================================
import os, re, copy, time, pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, arviz as az, itertools as it, xarray as xr
from copy import deepcopy
import matplotlib.collections as mpcollections
from scipy.stats import linregress
from arviz.utils import Numba
from Bio import Phylo, MissingPythonDependencyError
from io import StringIO
import argparse
# matplotlib.use('Agg')

def GetAllFiles(rootdir):
    filepath=[]
    sublist = os.listdir(rootdir)
    for sub in sublist:
        subdir = os.path.join(rootdir,sub)
        if os.path.isdir(subdir):
            filepath.extend(GetAllFiles(subdir))
        else:
            filepath.append(subdir)
    return filepath

def get_trees(stdouttxt):
    edge_line0 = [i for i in range(len(stdouttxt)) if stdouttxt[i].startswith(' father   node  name')][0] + 1
    edge_line1 = [i for i in range(len(stdouttxt)) if stdouttxt[i] == '\n' and i >= edge_line0][0]
    tree_edge = stdouttxt[edge_line0:edge_line1]
    tree_edge = [re.split(' +',i.replace('\n','').strip())[:2] for i in tree_edge]
    tree_edge = pd.DataFrame(tree_edge, columns=['parent', 'daughter'])
    number_tree = Phylo.read(StringIO(stdouttxt[edge_line1+1].replace('\n','')), 'newick', rooted=True)
    label_tree = Phylo.read(StringIO(stdouttxt[edge_line1+2].replace('\n','')), 'newick', rooted=True)
    for n in range(len(number_tree.get_nonterminals())):
        label_tree.get_nonterminals()[n].name = str(n+len(label_tree.get_terminals())+1)
        number_tree.get_nonterminals()[n].name = str(n+len(number_tree.get_terminals())+1)
    return label_tree, number_tree, tree_edge

def get_timepriors(stdouttxt):
    calibr_pattern = re.compile(r'^Node\s\d+:\s+.+,.+,.+,.+\)\n$')
    tmpr = [i.replace(' ','').replace('Node','').replace('\n','').replace(')','') for i in list(filter(calibr_pattern.search, stdouttxt))]
    tmpr = [re.split(r'[:\(,]', i)[:6] for i in tmpr]
    tmpr = pd.DataFrame(tmpr, columns = ['node', 'distribution', 'p1', 'p2', 'p3', 'p4'])
    tmpr.drop_duplicates(subset=['node'], inplace=True)
    tmpr.index = tmpr.iloc[:,0]
    tmpr = tmpr.iloc[:,1:]
    # tmpr.loc[:,'min'] = tmpr.loc[:,'min'].astype(float) * 100
    # tmpr.loc[:,'max'] = tmpr.loc[:,'max'].astype(float) * 100
    # tmpr.loc[tmpr['range']=='L', 'max'] = None
    # tmpr.loc[tmpr['range']=='U', 'max'] = tmpr.loc[tmpr['range']=='U', 'min']
    # tmpr.loc[tmpr['range']=='U', 'min'] = None
    return tmpr

def get_info_from_stdout(stdoutfilepath):
    with open(stdoutfilepath, 'r') as f:
        stdouttxt = f.readlines()
        label_tree, number_tree, tree_edge = get_trees(stdouttxt)
        timepriors = get_timepriors(stdouttxt)
    return label_tree, number_tree, tree_edge, timepriors

def create_inference_data_from_mcmc(postdir, priodir):
    # Assemble posterior and prior MCMC samples into an ArviZ InferenceData object.
    # Recursively search each directory for `mcmc.txt` files, read each chain,
    # add explicit chain and draw coordinates, and then store the data as
    # posterior and prior datasets for diagnostics and plotting.
    filepath = GetAllFiles(postdir)
    post_mcmcfilepath = [i.replace('\\','/') for i in filepath if i.endswith('mcmc.txt')]
    post_mcmcchains = list(map(lambda x: pd.read_csv(x, sep='\t', header=0, index_col=0), post_mcmcfilepath))
    chain = [i for i in range(len(post_mcmcchains)) for j in range(len(post_mcmcchains[i]))]
    draw = [j for i in range(len(post_mcmcchains)) for j in range(len(post_mcmcchains[i]))]
    post_mcmcchains = pd.concat(post_mcmcchains)
    post_mcmcchains['chain'] = chain
    post_mcmcchains['draw'] = draw
    post_mcmcchains.set_index(['chain', 'draw'], append=False, inplace=True)
    # get prior MCMC sample information from priodir mcmc files
    filepath = GetAllFiles(priodir)
    prio_mcmcfilepath = [i.replace('\\','/') for i in filepath if i.endswith('mcmc.txt')]
    prio_mcmcchains = list(map(lambda x: pd.read_csv(x, sep='\t', header=0, index_col=0), prio_mcmcfilepath))
    chain = [i for i in range(len(prio_mcmcchains)) for j in range(len(prio_mcmcchains[i]))]
    draw = [j for i in range(len(prio_mcmcchains)) for j in range(len(prio_mcmcchains[i]))]
    prio_mcmcchains = pd.concat(prio_mcmcchains)
    prio_mcmcchains['chain'] = chain
    prio_mcmcchains['draw'] = draw
    prio_mcmcchains.set_index(['chain', 'draw'], append=False, inplace=True)
    # combine posterior and prior MCMC samples into arviz-inference-data
    mcmcchains = az.InferenceData(posterior=xr.Dataset.from_dataframe(post_mcmcchains), prior=xr.Dataset.from_dataframe(prio_mcmcchains))
    return mcmcchains

def get_best_converged_chains(mcchains, r_hat_range=0.1, ess_min=100, mcse_max=0.1, all_chain_privileged=True, age_kind='median', at_least=2):
    # Select the best subset of chains based on convergence diagnostics.
    # Returns all chains when they all satisfy R-hat, ESS, and MCSE thresholds.
    # Otherwise, examine chain combinations and choose the best-performing subset.
    nchains = mcchains.sizes['chain']
    mcse_m = 'mcse_' + age_kind.lower()
    ess_m = 'ess_' + age_kind.lower()
    ess_m_r = ess_m + '_rel'
    def mcmcdiag4(x):
        y = az.summary(mcchains.sel(chain=x), kind='diagnostics', round_to='none', stat_focus=age_kind.lower())
        return pd.Series({'chain_num':x, 'nchs':len(x), mcse_m:max(y[mcse_m]), ess_m_r:min(y[ess_m]/len(x)), 'ess_tail_rel':min(y['ess_tail']/len(x)), 'r_hat_diff':max(y['r_hat'] - 1)})
    x4 = mcmcdiag4(list(range(nchains)))
    if all_chain_privileged and 0 <= x4['r_hat_diff'] < r_hat_range and x4[ess_m_r] >= ess_min and x4['ess_tail_rel'] >= ess_min and x4[mcse_m] <= mcse_max:
        return [ch for ch in range(nchains)], []
    else:
        mcmcdiag_all = pd.DataFrame({'chain_num':sum([[list(ch) for ch in it.combinations(range(nchains),i)] for i in range(at_least, nchains)],[])})
        mcmcdiag_all = mcmcdiag_all['chain_num'].apply(mcmcdiag4)
        mcmcdiag_all = pd.concat([mcmcdiag_all, x4.to_frame().T], axis=0, ignore_index=True)
        mcmcdiag_all.dropna(axis=0, how='any', subset=['r_hat_diff',ess_m_r,'ess_tail_rel',mcse_m], inplace=True)
        redline = pd.DataFrame([[0 <= rh < r_hat_range for rh in mcmcdiag_all['r_hat_diff']], [ess_min <= ess for ess in mcmcdiag_all[ess_m_r]], [ess_min <= ess for ess in mcmcdiag_all['ess_tail_rel']], [mcse < mcse_max for mcse in mcmcdiag_all[mcse_m]]])
        redline = redline.apply(all, axis=0)
        if redline.sum():
            mcmcdiag_all = mcmcdiag_all.loc[redline]
            mcmcdiag_all.sort_values(by=[ess_m_r, 'ess_tail_rel', 'nchs'], ignore_index=True, ascending=False, inplace=True)
        else:
            mcmcdiag_all.sort_values(by=['r_hat_diff', mcse_m], ignore_index=True, ascending=True, inplace=True)
        # mcmcdiag_all.drop_duplicates('nchs', keep='last', ignore_index=True, inplace=True)
        best_index = mcmcdiag_all['chain_num'][0]
        return [ch for ch in best_index], [ch for ch in range(nchains) if ch not in best_index]

def plot_mcmc_trace(mcmc, chains_used, filename_prefix, dpi=300):
    chaincolors = [
        "#0072B2FF",  # Blue
        "#E69F00FF",  # Orange
        "#009E73FF",  # Green
        "#CC79A7FF",  # Purple
        "#D55E00FF",  # Red
        "#F0D942FF",  # Yellow
        "#56B4E9FF",  # Cyan
        "#A6761DFF",  # Brown
        "#FF6347FF",  # Tomato
        "#4682B4FF"   # Steel Blue
    ]
    # # Create a figure that show color palette in preview
    # plt.figure(figsize=(10, 2))
    # for i, color in enumerate(chaincolors):
    #     plt.fill_between([i, i + 1], 0, 1, color=color, label=f"Color {i + 1}")
    # plt.xticks(range(len(chaincolors)), [f"{i + 1}" for i in range(len(chaincolors))])
    # plt.yticks([])
    # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    # plt.title("Nature/Science Color Palette (Colorblind-Friendly)")
    # plt.show()
    maxhight = len(mcmc.posterior.data_vars) * 2
    az.rcParams["plot.max_subplots"] = maxhight # set the maximum number of subplots larger than twice the number of parameters
    def plot_trace_group(group, chains, suffix):
        axes = az.plot_trace(
            group,
            figsize=(16, maxhight),
            compact=True,
            combined=False,
            chain_prop={"color": [chaincolors[i] for i in chains]},
            plot_kwargs={"alpha": 1, "linestyle": "-"},
            legend=True
        )
        # Adjust subplot layout
        plt.subplots_adjust(left=0.05, right=0.95, top=0.998, bottom=0.002, hspace=0.5, wspace=0.1)
        # Save the plot
        plt.savefig(f"{filename_prefix}{suffix}.pdf", dpi=dpi)
    # Plot posterior and prior groups
    plot_trace_group(mcmc.posterior, chains_used["posterior"], "mcmc_trace_posterior")
    plot_trace_group(mcmc.prior, chains_used["prior"], "mcmc_trace_prior")

def DateTree(label_tree, number_tree, tree_edge, node_info, calinodes, age_kind='median'):
    # Build a dated phylogenetic tree with branch lengths and node comments.
    # Uses estimated node ages to compute each branch length and attaches
    # calibration metadata for downstream Nexus export and plotting.
    NTips = label_tree.count_terminals()
    NNodes = len(label_tree.get_nonterminals())
    age = pd.concat([pd.Series([node_info[age_kind.lower()].max()] + [0.] * NTips, index=[str(i) for i in range(NTips+1)], name=age_kind.lower()), node_info[age_kind.lower()]])
    if age_kind.lower() == 'mean':
        cmt = lambda x: '&MEDIAN=' + str(x.iloc[0]) + ',MAD_RANGE={' + str(x.iloc[0]-x.iloc[1]) + ', ' + str(x.iloc[0]+x.iloc[1]) + '},95%ETI={' + str(x.iloc[2]) + ', ' + str(x.iloc[3]) + '},SD_RANGE={' + str(x.iloc[4]-x.iloc[5]) + ', ' + str(x.iloc[4]+x.iloc[5]) + '},95%HPD={' + str(x.iloc[6]) + ', ' + str(x.iloc[7]) + '},MCSE_MEAN=' + str(x.iloc[9]) + ',ESS_BULK=' + str(x.iloc[12]) + ',MCSE_MEDIAN=' + str(x.iloc[8]) + ',ESS_MEDIAN=' + str(x.iloc[11]) + ',ESS_TAIL=' + str(x.iloc[13]) + ',R_HAT=' + str(x.iloc[14])
    else:
        cmt = lambda x: '&MEAN=' + str(x.iloc[4]) + ',MAD_RANGE={' + str(x.iloc[0]-x.iloc[1]) + ', ' + str(x.iloc[0]+x.iloc[1]) + '},95%ETI={' + str(x.iloc[2]) + ', ' + str(x.iloc[3]) + '},SD_RANGE={' + str(x.iloc[4]-x.iloc[5]) + ', ' + str(x.iloc[4]+x.iloc[5]) + '},95%HPD={' + str(x.iloc[6]) + ', ' + str(x.iloc[7]) + '},MCSE_MEDIAN=' + str(x.iloc[8]) + ',ESS_MEDIAN=' + str(x.iloc[11]) + ',MCSE_MEAN=' + str(x.iloc[9]) + ',ESS_BULK=' + str(x.iloc[12]) + ',ESS_TAIL=' + str(x.iloc[13]) + ',R_HAT=' + str(x.iloc[14])
    node_cmt = node_info.apply(cmt, axis=1)
    branch_length = tree_edge.apply(lambda x: age[x.iloc[0]] - age[x.iloc[1]], axis=1)
    branch_length.index = tree_edge['daughter'].values
    time_tree = copy.deepcopy(label_tree)
    for nd in range(NNodes):
        time_tree.get_nonterminals()[nd].branch_length = branch_length.loc[number_tree.get_nonterminals()[nd].name]
        time_tree.get_nonterminals()[nd].comment = node_cmt.loc[number_tree.get_nonterminals()[nd].name]
        time_tree.get_nonterminals()[nd].name = time_tree.get_nonterminals()[nd].name + ['','@'][number_tree.get_nonterminals()[nd].name in calinodes] 
    for tp in range(NTips):
        time_tree.get_terminals()[tp].branch_length = branch_length.loc[number_tree.get_terminals()[tp].name]
    return time_tree

def TimePriPost(timepriors_test):
    # Compare prior and posterior calibration distributions for a single node.
    # This detects whether the posterior median is outside the prior HDI, and
    # whether the posterior 95% interval extends beyond the prior 95% bounds.
    if any(timepriors_test[['pr_hdi_97.5%','pr_hdi_2.5%']].isna()):
        return pd.Series([None, None])
    else:
        pos = ['OK', 'OK']
        Outflow = [timepriors_test['po_median'] > timepriors_test['pr_hdi_97.5%'], \
                timepriors_test['po_median'] < timepriors_test['pr_hdi_2.5%']]
        if Outflow[0]:
            pos[0] = 'U.out'
        elif Outflow[1]:
            pos[0] = 'L.out'
        else:
            pass    
        Outflow = [timepriors_test['po_hdi_97.5%'] < timepriors_test['pr_hdi_2.5%'], \
                   timepriors_test['po_hdi_2.5%'] > timepriors_test['pr_hdi_97.5%'], \
                   timepriors_test['po_hdi_97.5%'] > timepriors_test['pr_hdi_97.5%'], \
                   timepriors_test['po_hdi_2.5%'] < timepriors_test['pr_hdi_2.5%']]
        if Outflow[0]:
            pos[1] = 'L.out'
        elif Outflow[1]:
            pos[1] = 'U.out'
        elif all(Outflow[2:]):
            pos[1] = 'B.over'
        elif not any(Outflow[2:]):
            pass
        else:
            pos[1] = 'L.over' if Outflow[3] else 'U.over'
        return pd.Series(pos)

def iou(d1, d2, bins=1000):
    # Compute the Intersection over Union (IoU) between two sample distributions.
    # This is used to quantify overlap between posterior and prior marginal PDFs.
    d1 = np.array(d1)
    l1 = len(d1)
    d2 = np.array(d2)
    l2 = len(d2)
    d_range = (min(np.concatenate([d1,d2])), max(np.concatenate([d1,d2])))
    d1_bins, _ = np.histogram(d1, bins=bins, range=(d_range))
    d2_bins, _ = np.histogram(d2, bins=bins, range=(d_range))
    overlap = np.minimum(d1_bins / l1, d2_bins / l2)
    union = np.maximum(d1_bins / l1, d2_bins / l2)
    return overlap.sum() / union.sum()

def p_mcmc(mcmc):
    # Estimate posterior/prior overlap for each parameter using IoU.
    # Returns a vector of overlap scores that indicate how distinct the posterior
    # distribution is from the prior distribution for each node parameter.
    post_mcmc = mcmc.to_dataframe(groups='posterior', var_names=['~lnL']).iloc[:,2:]
    prio_mcmc = mcmc.to_dataframe(groups='prior').iloc[:,2:]
    if post_mcmc.shape[1] != prio_mcmc.shape[1]:
        raise ValueError('Check the number of parameters!')
    post_mcmc = np.transpose(np.array(post_mcmc))
    prio_mcmc = np.transpose(np.array(prio_mcmc))
    pmcmc = list(map(lambda x: iou(x[0], x[1]), zip(post_mcmc, prio_mcmc)))
    pmcmc.append(None)
    return np.array(pmcmc)

def ppifs_plot(data, filename):
    # Create the infinite-sites plot, comparing posterior HDI width with median age.
    # A strong correlation suggests sufficient data information, while a flat or
    # noisy trend may indicate limiting information for deep node age estimates.
    po_hdi_width = data['po_hdi_97.5%'] - data['po_hdi_2.5%']
    po_slope, po_intercept, po_r, po_p, po_std_err = linregress(data['po_median'], po_hdi_width)
    po_exp = po_slope * data['po_median'] + po_intercept # expression
    pr_hdi_width = data['pr_hdi_97.5%'] - data['pr_hdi_2.5%']
    pr_slope, pr_intercept, pr_r, pr_p, pr_std_err = linregress(data['pr_median'], pr_hdi_width)
    pr_exp = pr_slope * data['pr_median'] + pr_intercept # expression
    fig = plt.figure(figsize=(16,12))
    plt.scatter(data['po_median'], po_hdi_width, color='blue', marker='s', label='posterior')
    plt.scatter(data['pr_median'], pr_hdi_width, color='white', edgecolors='blue', label='prior')
    plt.plot(data['po_median'], po_exp, color="red", label='posterior', linewidth=3)
    plt.plot(data['pr_median'], pr_exp, color="red", linestyle='--', label='prior', linewidth=3)
    plt.title("Posterior and Prior Infinite-Sites Plot", fontsize=20)
    plt.xlabel('Median (100 My)', fontsize=18)
    plt.ylabel('HDI Width (100 My)', fontsize=18)
    plt.text(0.15, 0.85, f'Pos: y={round(po_slope, 3)} x+{round(po_intercept, 3)}, r2={round(po_r**2, 4)}, p={round(po_p, 4)}\nPri: y={round(pr_slope, 3)} x+{round(pr_intercept, 3)}, r2={round(pr_r**2, 4)}, p={round(pr_p, 4)}', ha='left', va='top', transform=fig.transFigure, fontsize=16)
    plt.legend()
    plt.savefig(f'{filename}.pdf', dpi=300)
    # plt.show()
    regstat = pd.DataFrame(data={'slope': [po_slope, pr_slope], 'std_err': [po_std_err, pr_std_err], 'intercept': [po_intercept, pr_intercept], 'r2': [po_r ** 2, pr_r ** 2], 'p_value': [po_p, pr_p]}, index=['posterior', 'prior'])
    return regstat

def pp_pdf_compare(mcmc, filename):
    # Plot prior and posterior marginal distributions for node ages.
    # This visual comparison helps assess how much the data updates each prior.
    nodelables = [i.replace('t_n','') for i in list(mcmc.posterior.data_vars) if i.startswith('t_n')]
    Nnodes = len(nodelables)
    post_mcmc = np.array(mcmc.to_dataframe(groups='posterior', var_names=['~lnL']).loc[:,[i for i in list(mcmc.posterior.data_vars) if i.startswith('t_n')]])
    prio_mcmc = np.array(mcmc.to_dataframe(groups='prior').loc[:,[i for i in list(mcmc.prior.data_vars) if i.startswith('t_n')]])
    Nsamples_post = post_mcmc.shape[0]
    Nsamples_prio = prio_mcmc.shape[0]
    post_mcmc = np.transpose(post_mcmc).reshape(-1,1)
    prio_mcmc = np.transpose(prio_mcmc).reshape(-1,1)
    comb_mcmc = pd.DataFrame(np.concatenate((post_mcmc, prio_mcmc)), columns=['mcmc'])
    comb_mcmc['Distribution'] = ['Posterior'] * (Nsamples_post * Nnodes) + ['Prior'] * (Nsamples_prio * Nnodes)
    comb_mcmc['node'] = np.concatenate((np.repeat(nodelables, repeats=Nsamples_post), np.repeat(nodelables, repeats=Nsamples_prio)))
    fig = plt.figure(figsize=(20,160))
    sns.violinplot(x='mcmc', y='node', data=comb_mcmc, hue='Distribution', split=True, dodge="auto", orient='h', linewidth = 1, width = 0.8, palette = 'Set2', order = nodelables, density_norm = 'count', gridsize = 30, gap=0, inner = 'quartile', bw_adjust = 0.01).invert_xaxis()
    plt.title("Posterior and Prior MCMC PDF Overlap", fontsize=24)
    plt.xlabel('Age (100 Ma)', fontsize=20)
    plt.ylabel('Node', fontsize=20)
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=20)
    plt.tight_layout()
    plt.savefig(f'{filename}.pdf', dpi=300)

def draw_timetree(tree, label_func=str, do_show=True, show_confidence=True, axes=None, branch_labels=None, node_bars=None, node_bar_color='blue', label_colors=None, font_size=None, *args, **kwargs,):
    # Draw a dated tree using matplotlib. Supports branch labels, optional
    # confidence values, custom label colors, and node bars for uncertainty.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MissingPythonDependencyError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None
    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []
    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)
    if not branch_labels:
        if show_confidence:
            def format_branch_label(clade):
                try:
                    confidences = clade.confidences
                    # phyloXML supports multiple confidences
                except AttributeError:
                    pass
                else:
                    return "/".join(conf2str(cnf.value) for cnf in confidences)
                if clade.confidence is not None:
                    return conf2str(clade.confidence)
                return None
        else:
            def format_branch_label(clade):
                return None
    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels    
    # optiond for node bars.
    if not node_bars:
        def format_node_bar(clade):
            return None
    elif type(node_bars) == str:
        def format_node_bar(clade):
            clade_cmt = eval('{'+clade.comment.replace('=','":').replace('&','"').replace(',',',"').replace('{','[').replace('}',']').replace('," ', ',')+'}')
            if node_bars in ['SD_RANGE','95%HPD']:
                clade_bar = clade_cmt.get(node_bars)
            return clade_bar
    elif callable(node_bars):
        def format_node_bar(clade):
            return node_bars(clade)
    # options for displaying label colors.
    if label_colors:
        if callable(label_colors):
            def get_label_color(label):
                return label_colors(label)
        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")
    else:
        def get_label_color(label):
            # if label_colors is not specified, use black
            return "black"
    # Layout
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.
        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths
    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.
        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
        }
        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0
        if tree.root.clades:
            calc_row(tree.root)
        return heights
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        if font_size:
            plt.rcparams["font.size"] = font_size["font.size"]
            plt.rcparams["axes.titlesize"] = font_size["axes.titlesize"]
            plt.rcparams["axes.labelsize"] = font_size["axes.labelsize"]
            plt.rcparams["xtick.labelsize"] = font_size["xtick.labelsize"]
            plt.rcparams["ytick.labelsize"] = font_size["ytick.labelsize"]
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError(f"Invalid argument for axes: {axes}")
    def draw_clade_lines(use_linecollection=False, orientation="horizontal", y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0, color="black", lw=".1",):
        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection([[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw))
        elif not use_linecollection and orientation == "vertical":
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection([[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw))
    def draw_clade(clade, x_start, color, lw, node_bars, node_bar_color,):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade] - max(x_posns.values())
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(use_linecollection=True, orientation="horizontal", y_here=y_here, x_start=x_start, x_here=x_here, color=color, lw=lw,)
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            axes.text(x_here, y_here, f" {label}", verticalalignment="center", color=get_label_color(label),)
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            axes.text(0.5 * (x_start + x_here), y_here, conf_label, fontsize="large", horizontalalignment="center",)
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation="vertical", x_here=x_here, y_bot=y_bot, y_top=y_top, color=color, lw=lw,)
            # Add node bars (optional)
            if node_bars:
                node_bar = format_node_bar(clade)
                axes.hlines(y_here, -node_bar[1], -node_bar[0], color=node_bar_color, lw=4*lw)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw, node_bars, node_bar_color)
    draw_clade(tree.root, -max(tree.depths().values()) * 1.01, "k", plt.rcParams["lines.linewidth"], node_bars, node_bar_color)
    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)
    # Aesthetics
    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            axes.set_title(name)
    axes.set_xlabel("Age (Ma)")
    axes.set_ylabel("")
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-1.05 * xmax, 1.2 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)
    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))
    if do_show:
        plt.show()

def set_colors(species_class):
    classes=species_class.drop_duplicates()
    colors=plt.cm.jet(np.linspace(0, 1, len(classes)))
    colors=pd.Series(dict(zip(classes, colors)))
    Xcolors=species_class.apply(lambda x: colors[x])
    return Xcolors, colors

def plot_phylo(tree, filename, width, height, tip_colorclass=None, bar_class=None,):
    # from phylogenetic_tree_draw import draw_timetree
    xtree = deepcopy(tree)
    max_age = max(xtree.depths().values())
    for i in xtree.get_nonterminals():
        i.name = ''
    fig = plt.figure(figsize=(width,height))
    axes = fig.add_subplot(1, 1, 1)
    fs = {"font.size":48, "axes.titlesize":52, "axes.labelsize":52, "xtick.labelsize":52, "ytick.labelsize":52, "legend.fontsize":36, "lines.linewidth":12}
    if isinstance(tip_colorclass, pd.DataFrame):
        draw_timetree(xtree, do_show=False, show_confidence=False, branch_labels=None, axes=axes, label_colors=tip_colorclass['Color'].to_dict(), node_bars='95%HPD', node_bar_color='blue', font_size = fs, )
        # tip color legend
        tipcolor_legend = tip_colorclass.drop_duplicates(inplace=False)
        tipcolor_legend = tipcolor_legend.sort_values(by='Color_label', ascending=True)
        for i in range(len(tipcolor_legend)):
            plt.text(x=-max_age*0.8, y=(xtree.count_terminals()-i*0.5), s=tipcolor_legend['Color_label'].iloc[i], fontsize=30, color=tipcolor_legend['Color'].iloc[i])
    else:
        draw_timetree(xtree, do_show=False, show_confidence=False, branch_labels=None, axes=axes, node_bars='95%HPD', node_bar_color='blue', font_size = fs, )
    # class bars
    if isinstance(bar_class, pd.DataFrame):
        bar_class = bar_class.loc[[i.name for i in xtree.get_terminals()]]
        Nclasses = bar_class.shape[1]
        for cl in range(Nclasses):
            lab = [bar_class.iloc[0,cl]]
            col = set_colors(bar_class.iloc[:,cl])[1]
            ypos = [[0.7, 1.3]]
            for i in range(1, len(bar_class.iloc[:,cl])):
                if bar_class.iloc[i,cl] == bar_class.iloc[i-1,cl]:
                    ypos[-1][1] = i + 1.3
                else:
                    lab.append(bar_class.iloc[i,cl])
                    ypos.append([i + 0.7, i + 1.3])
            for i in range(len(lab)):
                plt.text(x=max_age*0.4+max_age*0.8*(cl/Nclasses+0.01), y=ypos[i][0]+0.3, s=lab[i], va='center', fontsize=30, color=col[lab[i]])
                plt.plot([max_age*0.4+max_age*0.8*cl/Nclasses, max_age*0.4+max_age*0.8*cl/Nclasses], ypos[i], color=col[lab[i]], linewidth=6)
    # vertical scales per 100 my
    for i in range(round(max_age/100)+1):
        plt.axvline(x= -100 * i, color='grey', linewidth=2, linestyle='--')
    plt.title("Phylogenetic Time Tree", fontsize=28)
    plt.tight_layout()
    plt.savefig(f'{filename}.pdf', dpi=300)
    # plt.show()

if __name__ == "__main__":
    # Main entry point for the script when executed from the command line.
    # Parses directories, file names, and options, then runs the full MCMCtree
    # summary and visualization workflow.
    parser = argparse.ArgumentParser()
    parser.add_argument('--postdir', '-p', help='Directory containing the output files of posterior sampling.', required=True)
    parser.add_argument('--priodir', '-d', help='Directory containing the output files of prior sampling.', required=True)
    parser.add_argument('--outfile', '-o', nargs=2, help='Name of the std.out files. By default, the program will search the specified directory (including sub-directories) for the "step2_log.text" and "step3_log.text" files.', default=['step2_log.text', 'step3_log.text'], required=False)
    parser.add_argument('--mcmcfile', '-m', help='Name of the MCMC file. By default, the program will search the specified directory (including sub-directories) for the "mcmc.txt" files.', default='mcmc.txt', required=False)
    parser.add_argument('--classification', '-c', help='Name of the classification file. It should be an EXCEL file and is used to assign coloirs and vertical bars for grouping of tip labels.', default=None, required=False)
    parser.add_argument('--tipcolor', '-t', help='An excel file that specifying the color of tip labels, default tip color is "black".', default=None, required=False)
    parser.add_argument('--burnin', '-b', type=int, help='Disregard first N samples. Here is the real sample, NOT the setting in your original mcmctree control file but the LINES saved in the "mcmc.txt" file.', default=0, required=False)
    parser.add_argument('--resample', '-x', type=int, help='Resample every X samples (including the first after burnin) for all the chains. Default is NO resampling', default=1, required=False)
    parser.add_argument('--prefix', '-f', help='Add a prefix of "xxxx" to the output as "xxxx_filename". Default is NO prefix', default='', required=False)
    parser.add_argument('--select_chains', '-s', action='store_true', help='Ask the programe to select MCMC chains/runs as per r_hat value (1 < r_hat < 1.1), ess > N * good_ess, and mcse < 0.1.')
    parser.add_argument('--good_ess', '-e', type=float, help='The minimum of per-chain effective sample size for the good convergence of chains, default = 100 * number of chains. The acceptable (acpt) convergence defined by a half of this value.', default=100, required=False)
    parser.add_argument('--rm_mcmc', '-r', help='Ask the programe to remove the MCMC-files (mcmc.txt) and compress the rest files at the end to save disk space. "all" means deleting MCMC-files of all chains; "good", deleting MCMC-files if it is a good-tree result; "bad", deleting MCMC-files if it is a bad-tree result; "auto", deleting unconverged or unstationary and retaining converged and stationary chains for the posterior samples and deleting all MCMC-files for the prior samples if it is a good-tree result or deleting all MCMC-files if it is a bad-tree result; "none", keeping everything (default).', default='none')
    parser.add_argument('--compress', '-z', action='store_true', help='Ask the programe to compress all MCMCTree folders into three tar.gz-files.')
    args = parser.parse_args()

    # reading data
    curdir = os.getcwd().replace('\\','/')
    sys_delimiter = '\\' if os.name == 'nt' else '/'
    filename_prefix = curdir + sys_delimiter + args.prefix + '_' if args.prefix else ''

    if os.path.isabs(args.postdir):
        postdir = args.postdir
    else:
        postdir = curdir + '/' + args.postdir.replace('./','')

    if os.path.isabs(args.priodir):
        priodir = args.priodir
    else:
        priodir = curdir + '/' + args.priodir.replace('./','')   
    
    # get posterior MCMC sample information from postdir std.out files
    filepath = GetAllFiles(postdir)
    post_mcmcfilepath = [i.replace('\\','/') for i in filepath if i.endswith(args.mcmcfile)]
    post_stdoutfilepath = [i.replace('\\','/') for i in filepath if i.endswith(args.outfile[0])][0]
    label_tree, number_tree, tree_edge, timepriors = get_info_from_stdout(post_stdoutfilepath)
    post_mcmcchains = list(map(lambda x: pd.read_csv(x, sep='\t', header=0, index_col=0), post_mcmcfilepath))
    
    # get prior MCMC sample information from priodir std.out files
    filepath = GetAllFiles(priodir)
    prio_mcmcfilepath = [i.replace('\\','/') for i in filepath if i.endswith(args.mcmcfile)]
    prio_stdoutfilepath = [i.replace('\\','/') for i in filepath if i.endswith(args.outfile[1])][0]
    prio_mcmcchains = list(map(lambda x: pd.read_csv(x, sep='\t', header=0, index_col=0), prio_mcmcfilepath))

    # # prepare data for tuning (comment this block out before running the script)
    # postdir='G:/BristolProjects/Nematoda/codes/DIR_DIV'
    # priodir='G:/BristolProjects/Nematoda/codes/DIR_PRD'
    # filepath = GetAllFiles(postdir)
    # post_stdoutfilepath = [i.replace('\\','/') for i in filepath if i.endswith('step2_log.text')][0]
    # label_tree, number_tree, tree_edge, timepriors = get_info_from_stdout(post_stdoutfilepath)
    # filepath = GetAllFiles(priodir)
    # prio_stdoutfilepath = [i.replace('\\','/') for i in filepath if i.endswith('step3_log.text')][0]

    # create InferenceData object from MCMC samples
    mcmcchains = create_inference_data_from_mcmc(postdir, priodir)

    # get posterior MCMC sample information from postdir std.out files 
    Ori_chains = mcmcchains.posterior.sizes['chain']
    Ori_samples = mcmcchains.posterior.sizes['draw']
    Ori_parameters = len(mcmcchains.posterior.data_vars) - 1
    Ori_treesize = label_tree.count_terminals()
    Ori_calibrations = len(timepriors)

    # burnin and reampling (not necessary as MCMCtree save samples after primary burnin and resampling)
    mcmcchains.sel(draw=slice(args.burnin, Ori_samples, args.resample), inplace=True)
    
    # find best converged samples
    if args.select_chains:
        post_runs_used, post_bad_runs = get_best_converged_chains(mcmcchains.posterior, r_hat_range=0.1, ess_min=args.good_ess * 0.5, mcse_max=0.1, all_chain_privileged=False, age_kind='median', at_least=2)
        mcmcchains.posterior = mcmcchains.posterior.sel(chain=post_runs_used)
        prio_runs_used, prio_bad_runs = get_best_converged_chains(mcmcchains.prior, r_hat_range=0.1, ess_min=args.good_ess * 0.5, mcse_max=0.1, all_chain_privileged=False, age_kind='median', at_least=2)
        mcmcchains.prior = mcmcchains.prior.sel(chain=prio_runs_used)
        
    # Final MCMC info
    Final_chains_post = mcmcchains.posterior.sizes['chain']
    Final_samples_post = mcmcchains.posterior.sizes['draw']
    Final_parameters_post = len(mcmcchains.posterior.data_vars) - 1
    Final_retained_chains_post = [x + 1 for x in post_runs_used]
    Pooled_samples_post = Final_samples_post * Final_chains_post

    Final_chains_prio = mcmcchains.prior.sizes['chain']
    Final_samples_prio = mcmcchains.prior.sizes['draw']
    Final_parameters_prio = len(mcmcchains.prior.data_vars) - 1
    Final_retained_chains_prio = [x + 1 for x in prio_runs_used]
    Pooled_samples_prio = Final_samples_prio * Final_chains_prio

    # Plot MCMC trace
    plot_mcmc_trace(mcmc=mcmcchains, chains_used={'posterior':post_runs_used, 'prior':prio_runs_used}, filename_prefix=filename_prefix, dpi=300)

    # MCMC summary
    post_SM = az.summary(mcmcchains.posterior, hdi_prob=0.95, round_to='none', stat_focus='median')
    post_SM_1 = az.summary(mcmcchains.posterior, hdi_prob=0.95, round_to='none', stat_focus='mean')
    post_SM = pd.concat([post_SM.iloc[:, :4], post_SM_1.iloc[:, :4], post_SM.iloc[:, 4:5], post_SM_1.iloc[:, 4:6], post_SM.iloc[:, 5:6], post_SM_1.iloc[:, 6:]], axis=1)
    post_SM.iloc[:, 8:] = post_SM.iloc[:, 8:].round(4)
    post_SM.index = [i.replace('x[', '').replace(']','') for i in post_SM.index]
    prio_SM = az.summary(mcmcchains.prior, hdi_prob=0.95, round_to='none', stat_focus='median')
    prio_SM_1 = az.summary(mcmcchains.prior, hdi_prob=0.95, round_to='none', stat_focus='mean')
    prio_SM = pd.concat([prio_SM.iloc[:, :4], prio_SM_1.iloc[:, :4], prio_SM.iloc[:, 4:5], prio_SM_1.iloc[:, 4:6], prio_SM.iloc[:, 5:6], prio_SM_1.iloc[:, 6:]], axis=1)
    prio_SM.index = [i.replace('x[', '').replace(']','') for i in prio_SM.index]
    prio_SM.iloc[:, 8:] = prio_SM.iloc[:, 8:].round(4)
    
    # Re-Date the phylogeny according to the combined MCMC samples
    node_info = prio_SM[[i.startswith('t_n') for i in prio_SM.index]]
    node_info = pd.concat([node_info.iloc[:, :8] * 100, node_info.iloc[:, 8:]], axis=1)
    node_info.index = [i.replace('t_n','') for i in node_info.index]
    combined_prior_tree = DateTree(label_tree, number_tree, tree_edge, node_info, calinodes=timepriors.index, age_kind='median')
    node_info = post_SM[[i.startswith('t_n') for i in post_SM.index]]
    node_info = pd.concat([node_info.iloc[:, :8] * 100, node_info.iloc[:, 8:]], axis=1)
    node_info.index = [i.replace('t_n','') for i in node_info.index]
    combined_posterior_tree = DateTree(label_tree, number_tree, tree_edge, node_info, calinodes=timepriors.index, age_kind='median')
    post_SM.columns = 'po_' + post_SM.columns
    prio_SM.columns = 'pr_' + prio_SM.columns
    SM = pd.concat([post_SM, prio_SM], axis=1)

    # Test prior vs. posterior distributions
    SM[['median.pos','hdi.pos']] = SM.apply(TimePriPost, axis=1)
    # SM.loc['lnL',['po_ess_tail','po_r_hat']] = None
    SM['p.mcmc'] = p_mcmc(mcmc=mcmcchains)
    SM['calibrated_node'] = ['*' if any([x.endswith(str(y)) for y in timepriors.index]) else None for x in SM.index]

    # visualisation of the combined MCMC samples
    # Infinite-sites plot for posterior and prior clock dating
    ifsp_reg = ppifs_plot(SM[[i.startswith('t_n') for i in SM.index]], filename=filename_prefix + 'infinite_sites_plot')

    # Plot of the posterior and prior distributions
    pp_pdf_compare(mcmc=mcmcchains, filename=filename_prefix + 'posterior_prior_disributions')

    SM.to_excel(filename_prefix + 'mcmc_combined_summary.xlsx', index_label='parameter')
    node_info.to_excel(filename_prefix + 'node_age_ma_summary.xlsx', index_label='node')
    timepriors.to_excel(filename_prefix + 'time_prior_test.xlsx', index_label='node')
    ifsp_reg.to_excel(filename_prefix + 'infinite_sites_regression.xlsx', index_label='node')
    del post_SM, prio_SM, post_mcmcchains, prio_mcmcchains, node_info

    # Distinguish resulting time trees as per the convergence of MCMC chains
    # As per Gelman et al. (2014) Bayesian Data Analysis [p.287], threshold of r_hat is 1.1 and ESS should be greater than 10 * number of chains. Here I use higher standards for ESS.
    ESS_min_post = args.good_ess * Final_chains_post
    if all(list(1 <= rh < 1.1 for rh in SM[:-1]['po_r_hat']) + list(SM[:-1]['po_mcse_median'] < 0.1) + [SM[:-1]['po_ess_median'].min() >= ESS_min_post] + [SM[:-1]['po_ess_tail'].min() >= ESS_min_post]):
        treemark_po='good'    
    elif all(list(1 <= rh < 1.1 for rh in SM[:-1]['po_r_hat']) + list(SM[:-1]['po_mcse_median'] < 0.1) + [SM[:-1]['po_ess_median'].min() >= ESS_min_post * 0.5] + [SM[:-1]['po_ess_tail'].min() >= ESS_min_post * 0.5]):
        treemark_po='acpt'
    else:
        treemark_po='bad'

    ESS_min_prio = args.good_ess * Final_chains_prio
    if all(list(1 <= rh < 1.1 for rh in SM[:-1]['pr_r_hat']) + list(SM[:-1]['pr_mcse_median'] < 0.1) + [SM[:-1]['pr_ess_median'].min() >= ESS_min_prio] + [SM[:-1]['pr_ess_tail'].min() >= ESS_min_prio]):
        treemark_pr='good'    
    elif all(list(1 <= rh < 1.1 for rh in SM[:-1]['pr_r_hat']) + list(SM[:-1]['pr_mcse_median'] < 0.1) + [SM[:-1]['pr_ess_median'].min() >= ESS_min_prio * 0.5] + [SM[:-1]['pr_ess_tail'].min() >= ESS_min_prio * 0.5]):
        treemark_pr='acpt'
    else:
        treemark_pr='bad'

    # Save results
    Phylo.write(combined_posterior_tree, filename_prefix + 'combined_posterior_tree_' + treemark_po + '.tre', "nexus")
    with open(filename_prefix + 'combined_posterior_tree_' + treemark_po + '.tre', 'r') as f:
        unrootedtree = f.read()
        rootedtree = unrootedtree.replace('Begin Trees;\n Tree tree1=(', 'Begin Trees;\n Tree tree1=[&R]((').replace('];\nEnd;\n', ']);\nEnd;\n')
    with open(filename_prefix + 'combined_posterior_tree_' + treemark_po + '.tre', 'w') as f:
        f.write(rootedtree)

    Phylo.write(combined_prior_tree, filename_prefix + 'combined_prior_tree_' + treemark_pr + '.tre', "nexus")
    with open(filename_prefix + 'combined_prior_tree_' + treemark_pr + '.tre', 'r') as f:
        unrootedtree = f.read()
        rootedtree = unrootedtree.replace('Begin Trees;\n Tree tree1=(', 'Begin Trees;\n Tree tree1=[&R]((').replace('];\nEnd;\n', ']);\nEnd;\n')
    with open(filename_prefix + 'combined_prior_tree_' + treemark_pr + '.tre', 'w') as f:
        f.write(rootedtree)

    # plot classification trees
    if args.classification is not None:
        sp_class = pd.read_excel(args.classification if os.path.isabs(args.classification) else curdir + '/' + args.classification.replace('./',''), header=0, index_col=0).sort_index()
    else:
        sp_class = None
    if args.tipcolor is not None:
        tip_colors = pd.read_excel(args.tipcolor if os.path.isabs(args.tipcolor) else curdir + '/' + args.tipcolor.replace('./',''), header=0, index_col=0).sort_index()
    else:
        tip_colors = None
       
    plot_phylo(combined_posterior_tree, filename=filename_prefix + 'combined_posterior_tree_' + treemark_po, width=60, height=160, tip_colorclass=tip_colors, bar_class=sp_class)

    plot_phylo(combined_prior_tree, filename=filename_prefix + 'combined_prior_tree_' + treemark_pr, width=60, height=160, tip_colorclass=tip_colors, bar_class=sp_class)    
    
    # Report
    report = 'MCMCtree Final Report\n' \
        + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) \
        + '\n--------------------\n' \
        + 'MCMCtree results summary\n' \
        + 'N. chains: ' + str(Ori_chains) + '\n' \
        + 'N. samples per chain: ' + str(Ori_samples) + '\n' \
        + 'N. parameters (excl. Ln-likelihood): ' + str(Ori_parameters) + '\n' \
        + 'N. calibrated nodes: ' + str(Ori_calibrations) + '\n' \
        + 'Tree size: ' + str(Ori_treesize) + '\n' \
        + '\n--------------------\n' \
        + 'Resampling\n' \
        + 'N. burnin per chain: ' + str(args.burnin) + '\n' \
        + 'Resampling interval: ' + str(args.resample) + '\n' \
        + '\n--------------------\n' \
        + 'Final retained MCMC samples\n' \
        + 'N. chains (posterior): ' + str(Final_chains_post) + '\n' \
        + 'N. samples per chain (posterior): ' + str(Final_samples_post) + '\n' \
        + 'N. retained chains (posterior): ' + str(Final_retained_chains_post) + '\n' \
        + 'N. chains (prior): ' + str(Final_chains_prio) + '\n' \
        + 'N. samples per chain (prior): ' + str(Final_samples_prio) + '\n' \
        + 'N. retained chains (prior): ' + str(Final_retained_chains_prio) + '\n' \
        + 'N. parameters (excl. Ln-likelihood): ' + str(Final_parameters_post) + '\n' \
        + 'N. calibrated nodes: ' + str(Ori_calibrations) + '\n' \
        + 'Tree size: ' + str(Ori_treesize) + '\n' \
        + '\n--------------------\n' \
        + 'Combined MCMC statistics (posterior)\n' \
        + 'N. pooled samples: ' + str(Pooled_samples_post) + '\n' \
        + 'Log Likelihood (lnL-Median[med. abs. deviation]): ' + str(SM.loc['lnL']['po_median']) + ' [' + str(SM.loc['lnL']['po_mad'].round(4)) + ']\n' \
        + 'Log Likelihood (lnL-HDI95%): [' + str(SM.loc['lnL']['po_hdi_2.5%']) + ', ' + str(SM.loc['lnL']['po_hdi_97.5%']) + ']\n' \
        + 'R^ range (converged if close to 1, here [1, 1.1)): ' + str(SM[:-1]['po_r_hat'].min()) + ', ' + str(SM[:-1]['po_r_hat'].max()) + '\n' \
        + 'N. parameters with R^ >= 1.1 or < 1: ' + str(sum((SM[:-1]['po_r_hat'] >= 1.1) | (SM[:-1]['po_r_hat'] < 1))) + '\n' \
        + 'Max. MCSE (converged if close to 0, here < 0.1): ' + str(SM[:-1]['po_mcse_mean'].max()) + '\n' \
        + 'N. parameters with MCSE >= 0.1: ' + str(sum(SM[:-1]['po_mcse_mean'] >= 0.1)) + '\n' \
        + 'Min. Median Effective Sample Size (ESS_Median): ' + str(SM[:-1]['po_ess_median'].min()) + '\n' \
        + 'N. parameters with ESS_Median < ' + str(ESS_min_post * 0.5) + ' (acceptable >= ' + str(args.good_ess * 0.5) + ' * N. chains): ' + str(sum(SM[:-1]['po_ess_median'] < ESS_min_post * 0.5)) + '\n' \
        + 'N. parameters with ESS_Median < ' + str(ESS_min_post) + ' (good >= ' + str(args.good_ess) + ' * N. chains): ' + str(sum(SM[:-1]['po_ess_median'] < ESS_min_post)) + '\n' \
        + 'Min. Tail Effective Sample Size (ESS_Tail): ' + str(SM['po_ess_tail'].min()) + '\n' \
        + 'N. parameters with ESS_Tail < ' + str(ESS_min_post * 0.5) + ' (acceptable >= ' + str(args.good_ess * 0.5) + ' * N. chains): ' + str(sum(SM[:-1]['po_ess_tail'] < ESS_min_post * 0.5)) + '\n' \
        + 'N. parameters with ESS_Tail < ' + str(ESS_min_post) + ' (good >= ' + str(args.good_ess) + ' * N. chains): ' + str(sum(SM[:-1]['po_ess_tail'] < ESS_min_post)) + '\n' \
        + '\n--------------------\n' \
        + 'Combined MCMC statistics (prior)\n' \
        + 'N. pooled samples: ' + str(Pooled_samples_prio) + '\n' \
        + 'R^ range (converged if close to 1, here [1, 1.1)): ' + str(SM[:-1]['pr_r_hat'].min()) + ', ' + str(SM[:-1]['pr_r_hat'].max()) + '\n' \
        + 'N. parameters with R^ >= 1.1 or < 1: ' + str(sum((SM[:-1]['pr_r_hat'] >= 1.1) | (SM[:-1]['pr_r_hat'] < 1))) + '\n' \
        + 'Max. MCSE (converged if close to 0, here < 0.1): ' + str(SM[:-1]['pr_mcse_mean'].max()) + '\n' \
        + 'N. parameters with MCSE >= 0.1: ' + str(sum(SM[:-1]['pr_mcse_mean'] >= 0.1)) + '\n' \
        + 'Min. Median Effective Sample Size (ESS_Median): ' + str(SM[:-1]['pr_ess_median'].min()) + '\n' \
        + 'N. parameters with ESS_Median < ' + str(ESS_min_prio * 0.5) + ' (acceptable >= ' + str(args.good_ess * 0.5) + ' * N. chains): ' + str(sum(SM[:-1]['pr_ess_median'] < ESS_min_prio * 0.5)) + '\n' \
        + 'N. parameters with ESS_Median < ' + str(ESS_min_prio) + ' (good >= ' + str(args.good_ess) + ' * N. chains): ' + str(sum(SM[:-1]['pr_ess_median'] < ESS_min_prio)) + '\n' \
        + 'Min. Tail Effective Sample Size (ESS_Tail): ' + str(SM[:-1]['pr_ess_tail'].min()) + '\n' \
        + 'N. parameters with ESS_Tail < ' + str(ESS_min_prio * 0.5) + ' (acceptable >= ' + str(args.good_ess * 0.5) + ' * N. chains): ' + str(sum(SM[:-1]['pr_ess_tail'] < ESS_min_prio * 0.5)) + '\n' \
        + 'N. parameters with ESS_Tail < ' + str(ESS_min_prio) + ' (good >= ' + str(args.good_ess) + ' * N. chains): ' + str(sum(SM[:-1]['pr_ess_tail'] < ESS_min_prio)) + '\n' \
        + '\n--------------------\n' \
        + 'Time prior test\n' \
        + 'N. time priors (B, U, L, ST, G): ' + str(Ori_calibrations) + ' (' + str(sum(timepriors['distribution'] == 'B')) + ', ' + str(sum(timepriors['distribution'] == 'U')) + ', ' + str(sum(timepriors['distribution'] == 'L')) + ', ' + str(sum(timepriors['distribution'] == 'ST')) + ', ' + str(sum(timepriors['distribution'] == 'G')) + ')\n' \
        + 'N. calibrated nodes with posterior mean outside prior 95%HDI: ' + str(sum(SM.loc[['t_n'+str(i) for i in timepriors.index]]['median.pos'] != 'OK')) + '\n' \
        + 'N. calibrated nodes with 95%HDI over prior 95%HDI bounds: ' + str(sum(SM.loc[['t_n'+str(i) for i in timepriors.index]]['hdi.pos'] != 'OK')) + '\n' \
        + 'N. calibrated nodes with p.mcmc < 0.05 (posterior and prior overlap): ' + str(sum(SM.loc[['t_n'+str(i) for i in timepriors.index]]['p.mcmc'] < 0.05)) + '\n' \
        + '\n--------------------\n' \
        + 'Posterior vs. prior parameters\n' \
        + 'N. parameters with posterior mean outside prior 95%HDI: ' + str(sum(SM[:-1]['median.pos'] != 'OK')) + '\n' \
        + 'N. parameters with 95%HDI over prior 95%HDI bounds: ' + str(sum(SM[:-1]['hdi.pos'] != 'OK')) + '\n' \
        + 'N. parameters with p.mcmc < 0.05 (posterior and prior overlap): ' + str(sum(SM[:-1]['p.mcmc'] < 0.05)) + '\n' \
        + 'N. nodes with posterior mean age outside prior 95%HDI: ' + str(sum(SM.loc[[i.startswith('t_n') for i in SM.index]]['median.pos'] != 'OK')) + '\n' \
        + 'N. nodes with age 95%HDI over prior 95%HDI bounds: ' + str(sum(SM.loc[[i.startswith('t_n') for i in SM.index]]['hdi.pos'] != 'OK')) + '\n' \
        + 'N. nodes with p.mcmc < 0.05 (posterior and prior overlap): ' + str(sum(SM.loc[[i.startswith('t_n') for i in SM.index]]['p.mcmc'] < 0.05)) + '\n' \
        + '\n--------------------\n' \
        + 'HDI-width vs. posterior mean regression (Infinite-sites plot)\n' \
        + str(ifsp_reg) + '\n' \
        + '\n--------------------\n' \
        + 'Detail files saved in:\n' \
        + filename_prefix + 'combined_posterior_tree_' + treemark_po + '.tre\n' \
        + filename_prefix + 'combined_posterior_tree_' + treemark_po + '.pdf\n' \
        + filename_prefix + 'combined_prior_tree_' + treemark_pr + '.tre\n' \
        + filename_prefix + 'combined_prior_tree_' + treemark_pr + '.pdf\n' \
        + filename_prefix + 'mcmc_combined_summary.xlsx\n' \
        + filename_prefix + 'node_age_ma_summary.xlsx\n' \
        + filename_prefix + 'time_prior_test.xlsx\n' \
        + filename_prefix + 'infinite_sites_regression.xlsx\n' \
        + filename_prefix + 'infinite_sites_plot.pdf\n' \
        + filename_prefix + 'posterior_prior_disributions.pdf\n' \
        + filename_prefix + 'mcmc_trace_posterior.pdf\n' \
        + filename_prefix + 'mcmc_trace_prior.pdf\n' \
        + '\n--------End---------\n'
    with open(filename_prefix + 'MCMCtree_final_report.txt', 'w') as f:
        f.write(report) 

    # remove MCMC-files
    if args.rm_mcmc.lower == 'all':
        [os.remove(i) for i in post_mcmcfilepath]
        [os.remove(i) for i in prio_mcmcfilepath]
    elif args.rm_mcmc.lower == 'good':
        if treemark_po != 'bad':
            [os.remove(i) for j in post_runs_used for i in post_mcmcfilepath if i.endswith('run' + str(j) + '/mcmc.txt')]
        if treemark_pr != 'bad':
            [os.remove(i) for j in prio_runs_used for i in prio_mcmcfilepath if i.endswith('run' + str(j) + '/mcmc.txt')]
    elif args.rm_mcmc.lower == 'bad':
        [os.remove(i) for j in post_bad_runs for i in post_mcmcfilepath if i.endswith('run' + str(j) + '/mcmc.txt')]
        [os.remove(i) for j in prio_bad_runs for i in prio_mcmcfilepath if i.endswith('run' + str(j) + '/mcmc.txt')]
    elif args.rm_mcmc.lower == 'auto':
        [os.remove(i) for j in post_bad_runs for i in post_mcmcfilepath if i.endswith('run' + str(j) + '/mcmc.txt')]
        [os.remove(i) for i in prio_mcmcfilepath]
    else:
        pass
    
    # compress the rest files into three zip-files (or tar.gz-files)
    # I found in some cases tar-compression may cause issues, so here I use zip instead of tar.gz for better compatibility across different operating systems, but you can switch to tar.gz if you prefer.
    if args.compress:
        os.system('zip -r DIR_BV.zip DIR_BV && rm -r DIR_BV')
        os.system('zip -r DIR_DIV.zip DIR_DIV && rm -r DIR_DIV')
        os.system('zip -r DIR_PRD.zip DIR_PRD && rm -r DIR_PRD')
        # os.system('tar -zcvf DIR_BV.tar.gz DIR_BV && rm -r DIR_BV')
        # os.system('tar -zcvf DIR_DIV.tar.gz DIR_DIV && rm -r DIR_DIV')
        # os.system('tar -zcvf DIR_PRD.tar.gz DIR_PRD && rm -r DIR_PRD')
    