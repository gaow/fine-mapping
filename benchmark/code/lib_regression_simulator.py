import numpy as np
import os, copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pprint import pformat
from collections import OrderedDict

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError(item)

    def __deepcopy__(self, memo):
        return dotdict(copy.deepcopy(dict(self)))
    
class RegressionData(dotdict):
    def __init__(self, X = None, Y = None, Z = None):
        # FIXME: check if inputs are indeed numpy arrays
        self.debug = dotdict()
        self.x_mean = self.y_mean = self.z_mean = None
        self.X = X
        self.Y = Y
        self.Z = Z
        self.xcorr = None

    def get_summary_stats(self):
        '''
        Computer univariate regression for every X_j (N by 1) and Y_r (N by 1)
        Bhat: J by R matrix of estimated effects
        Shat: J by R matrix of SE of Bhat
        '''
        if self.Z is not None:
            self.remove_covariates()
        # Compute betahat
        XtX_vec = np.einsum('ji,ji->i', self.X, self.X)
        self.Bhat = (self.X.T @ self.Y) / XtX_vec[:,np.newaxis]
        # Compute se(betahat)
        Xr = self.Y - np.einsum('ij,jk->jik', self.X, self.B)
        Re = np.einsum('ijk,ijk->ik', Xr, Xr)
        self.Shat = np.sqrt(Re / XtX_vec[:,np.newaxis] / (self.X.shape[0] - 2))

    def remove_covariates(self):
        if self.Z is not None:
            self.Y -= self.Z @ (np.linalg.inv(self.Z.T @ self.Z) @ self.Z.T @ self.Y)
            self.Z = None

    def center_data(self):
        # for np.array: np.mean(Z, axis=0, keepdims=True)
        # for np.matrix, no keepdims argument
        if self.X is not None and self.x_mean is None:
            self.x_mean = np.mean(self.X, axis=0)
            self.X -= self.x_mean
        if self.Y is not None and self.y_mean is None:
            self.y_mean = np.mean(self.Y, axis=0)
            self.Y -= self.y_mean
        if self.Z is not None and self.z_mean is None:
            self.z_mean = np.mean(self.Z, axis=0)
            self.Z -= self.z_mean

    def set_xcorr(self, xcorr=None):
        if xcorr is not None:
            self.xcorr = np.array(xcorr)
        else:
            self.xcorr = np.corrcoef(self.X, rowvar = False)
            self.xcorr = (np.square(self.xcorr) * np.sign(self.xcorr)).astype(np.float16)

    def plot_xcorr(self, out, limit = 5000, size = 15):
        if isinstance(limit, tuple):
            start = max(0, limit[0])
            end = min(self.xcorr.shape[0], limit[1])
            xcorr = self.xcorr[start:end,start:end]
        else:
            # the correlation matrix
            limit = min(self.xcorr.shape[0], limit)
            xcorr = self.xcorr[0:limit,0:limit]
        # Generate a mask for the upper triangle
        mask = np.zeros_like(xcorr, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True
        fig, ax = plt.subplots(figsize=(size,size))
        if out.endswith('pdf'):
            raise ValueError('Please use png extension for output file.')
        print(f'Plotting figure {out} for {limit} markers (default limit set to 5000) ...')
        use_abs = np.sum(xcorr < 0) == 0
        cmap = sns.diverging_palette(250, 15, as_cmap=True)
        sns.heatmap(xcorr, ax = ax, mask=mask, cmap = cmap, vmin=-1 if not use_abs else 0,
                    vmax=1, square=True, xticklabels = False, yticklabels = False, 
                    linewidths=.5, cbar_kws={"shrink": .5}, center=0)  
        ax = plt.gca()
        print(f'Saving figure {out} ...')
        plt.savefig(out, dpi = 500)
        
    def permute_X_columns(self):
        '''
        Permute X columns, i.e. break blocked correlation structure
        '''
        np.random.shuffle(self.X)
        
    def plot_property_vector(self, yaxis, zaxis, xz_cutoff = None, out = '/tmp/1.png',
                            conf = {'title': '', 'ylabel': '', 'zlabel': ''}):
        '''
        - yaxis can be eg $\beta$ or log10BF or -log10Prob
        - zaxis can be some other quantity whose value will be 
        reflected by color shade
        - xz_cutoff: (c1, c2). c1 is correlation cutoff to highlight
        when c2 is satisfied by a given position on x-axis
        '''
        xaxis = [x+1 for x in range(len(yaxis))]
        cmap = sns.cubehelix_palette(start=2.8, rot=.1, as_cmap=True)
        f, ax = plt.subplots(figsize=(18,5))
        if zaxis is not None:
            points = ax.scatter(xaxis, yaxis, c=zaxis, cmap=cmap)
            f.colorbar(points, label=conf['zlabel'])
        else:
            points = ax.scatter(xaxis, yaxis, cmap=cmap)
        if xz_cutoff is not None and zaxis is not None:
            c1, c2 = xz_cutoff
            if len([i for i in zaxis if i > c2]) > 100:
                print('Too many to highlight!')
            else:
                for idx, item in enumerate(zaxis):
                    if item > c2:
                        ax.scatter(xaxis[idx], yaxis[idx], s=80, 
                                   facecolors='none', edgecolors='r')
                        for ii, xx in enumerate(self.xcorr[idx,:]):
                            if xx > c1 and xx < 1.0:
                                ax.scatter(xaxis[ii], yaxis[ii], 
                                           color='y', marker='+')
        ax.set_title(conf['title'])
        ax.set_ylabel(conf['ylabel'])
        plt.gca()
        plt.savefig(out, dpi = 500)
        
    def get_representative_features(self, block_r2 = 0.8, block_size = 10, max_indep_r2 = 0.02):
        '''
        Based on xcorr matrix, select "most representative features". 
        That is, these features are potentially most convoluted by other features (have stronger xcorr)
        yet are independent among each other.
        - block_r2: definition of correlated block -- abs squared correlation have to be > cutoff1
        - block_size: define a large enough block -- block size have to be > block_size
        - max_indep_r2: now select features that are completely independent -- r2 < max_indep_r2
        '''
        if self.xcorr is None:
            self.set_xcorr(None)
        # get r2 summary
        r2 = pd.DataFrame(self.xcorr)
        strong_r2_count = ((np.absolute(r2) > block_r2) * r2).sum(axis = 0).sort_values(ascending = False)
        strong_r2_count = strong_r2_count[strong_r2_count > block_size]
        # filter by r2
        exclude = []
        for x in strong_r2_count.index:
            if x in exclude:
                continue
            for y in strong_r2_count.index:
                if y in exclude or y == x:
                    continue
                if np.absolute(r2[x][y]) > max_indep_r2:
                    exclude.append(y)
        return [x for x in strong_r2_count.index if not x in exclude]

    def __str__(self):
        return pformat(self.__dict__, indent = 4)
    
class ResidualCorrelation:
    def __init__(self, mode, dim = 1):
        self.mode = mode
        self.dim = dim
        
    def apply(self):
        if self.mode == 'identity':
            return self.set_identity()
        else:
            raise ValueError(f"Residual mode {self.mode} not implemented.")

    def set_identity(self):
        if self.dim > 1:
            # multivariate case
            return np.identity(self.dim)
        else:
            return 1
    
class UnivariateMixture:
    '''Simulated distributions of Stephens 2017 (ASH paper)'''
    def __init__(self, dim):
        self.size = dim
        self.pi0 = 0
        self.pis = []
        self.mus = []
        self.sigmas = []
        self.coef = []
        
    def set_vanilla(self, amplitude):
        self.pis = [1]
        self.mus = [0]
        self.sigmas = [amplitude]
        
    def set_pi0(self, pi0):
        self.pi0 = pi0
        
    def set_spiky(self):
        self.pis = [0.4,0.2,0.2,0.2]
        self.mus = [0,0,0,0]
        self.sigmas = [0.25,0.5,1,2]
    
    def set_near_normal(self):
        self.pis = [2/3,1/3]
        self.mus = [0,0]
        self.sigmas = [1,2]
        
    def set_flat_top(self):
        self.pis = [1/7] * 7
        self.mus = [-1.5, -1, -.5 , 0, .5, 1, 1.5]
        self.sigmas = [0.5] * 7
        
    def set_skew(self):
        self.pis = [1/4,1/4,1/3,1/6]
        self.mus = [-2,-1,0,1]
        self.sigmas = [2,1.5,1,1]
        
    def set_big_normal(self):
        self.pis = [1]
        self.mus = [0]
        self.sigmas = [4]

    def set_bimodal(self):
        self.pis = [0.5, 0.5]
        self.mus = [-2, 2]
        self.sigmas = [1, 1]
        
    def get_effects(self):
        '''
        beta ~ \pi_0\delta_0 + \sum \pi_i N(mu_i, sigma_i)
        '''
        sigmas = np.diag(self.sigmas)
        assert (len(self.pis), len(self.pis)) == sigmas.shape
        masks = np.random.multinomial(1, self.pis, size = self.size)
        mix = np.random.multivariate_normal(self.mus, sigmas, self.size)
        self.coef = np.sum(mix * masks, axis = 1) * np.random.binomial(1, 1 - self.pi0, self.size)
        
    def swap_top_effects(self, given_index):
        '''Set top effects to given indices
        One can specify index, or use the "top_index"
        generated by RegressionData.get_representative_features()
        '''
        given_index = np.array(given_index, dtype=int)
        nb = np.zeros(len(self.coef))
        beta = sorted(self.coef, key=abs, reverse=True)
        for idx in given_index:
            nb[idx] = beta.pop(0)
        np.random.shuffle(beta)
        for idx in range(len(nb)):
            if not idx in given_index:
                nb[idx] = beta.pop(0)
        assert len(beta) == 0
        self.coef = np.array(nb)
        
    def sparsify_effects(self, num_non_zero, top_only = False):
        '''
        top_only: only keep top `num_non_zero` effects
        '''
        if top_only:
            big_beta_index = [i[0] for i in sorted(enumerate(self.coef), key = lambda x: np.absolute(x[1]), reverse = True)]
            selected_index = big_beta_index[:min(len(big_beta_index), num_non_zero)]
        else:
            selected_index = np.random.choice(np.arange(self.size), size=num_non_zero)
        for j in range(self.size):
            if j not in selected_index:
                self.coef[j] = 0
                
    def get_y(self, regression_data, pve = None, sigma = None):
        if sigma is None and pve is None:
            raise ValueError('Need one of sigma or pve.')
        if not (pve > 0 and pve < 1):
            raise ValueError(f'PVE has to be between 0 and 1, not {pve}.')
        if pve is not None:
            genetic_var = np.var(np.dot(regression_data.X, self.coef.T))
            pheno_var = genetic_var / pve
            self.residual_variance = pheno_var - genetic_var
        y = np.dot(regression_data.X, self.coef.T) + np.random.normal(0, np.sqrt(self.residual_variance), regression_data.X.shape[0])
        # y.reshape(len(y), 1)
        return y.T
        
    def __str__(self):
        params = ' + '.join(["{} N({}, {}^2)".format(x,y,z) for x, y, z in zip(self.pis, self.mus, self.sigmas)])
        return '{:.3f} \delta_0 + {:.3f} [{}]'.format(self.pi0, 1 - self.pi0, params)
    
class MultivariateMixture:
    '''Implementation of multivariate regression based simulations (using mixture priors on effect size)'''
    def __init__(self, dim):
        self.J, self.R = dim
        self.pis = OrderedDict([('null', 0)])
        self.Us = OrderedDict([('null', np.zeros((self.R, self.R)))])
        self.mus = dict([('zeros', np.zeros(self.R))])
        self.coef = []
        self.grid = (0.2, 0.4, 0.6, 0.8)
        self._init_canonical()

    def set_pi0(self, pi0):
        self.pis['null'] = pi0
        
    def set_grid(self, grid):
        self.grid = grid
        
    def _init_canonical(self):
        '''
        U is a dict of 
        - "identity" for the identity (effects are independent among conditions);
        - "singletons" for the set of matrices with just one non-zero entry x_{jj} = 1 (j=1,...,R); (effect specific to condition j);
        - "equal_effects" for the matrix of all 1s (effects are equal among conditions);
        - "simple_het" for a set of matrices with 1s on the diagonal and all off-diagonal elements equal to pho; (effects are correlated among conditions).
        '''
        pho = [0.25, 0.5, 0.75]
        self.Us['identity'] = np.identity(self.R)
        for i in range(self.R):
            self.Us[f'singleton_{i+1}'] = np.diagflat([1 if idx == i else 0 for idx in range(self.R)])
        self.Us['equal_effects'] = np.ones((self.R, self.R))
        for idx, item in enumerate(sorted(pho)):
            self.Us[f'simple_het_{idx+1}'] = np.ones((self.R, self.R)) * item
            np.fill_diagonal(self.Us[f'simple_het_{idx+1}'], 1)

    def _set_single_component(self, u_id):
        self.pis[u_id] = 1 - self.pis['null']
        for k in self.Us:
            if k not in self.pis:
                self.pis[k] = 0

    def set_indep(self):
        '''
        All weights are on identity effects
        '''
        self._set_single_component('identity')
                
    def set_high_het(self):
        '''
        All weights are on simple_het_1
        '''
        self._set_single_component('simple_het_1')

    def set_mid_het(self):
        '''
        All weights are on simple_het_2
        '''
        self._set_single_component('simple_het_2')

    def set_low_het(self):
        '''
        All weights are on simple_het_3
        '''
        self._set_single_component('simple_het_3')

    def set_shared(self):
        '''
        All weights are on equal effects
        '''
        self._set_single_component('equal_effects')

    def set_singleton(self):
        '''
        Here I simulate a mixture of singleton matrices each condition having equal chance to be causal.
        '''
        weight = (1 - self.pis['null']) / self.R
        for item in range(self.R):
            self.pis[f'singleton_{item+1}'] = weight
        for k in self.Us:
            if k not in self.pis:
                self.pis[k] = 0

    def set_manual_mixture(self, cfg = dict(identity=0.1,equal_effects=0.2,singleton=0.2,simple_het_1=0.1,simple_het_2=0.1,simple_het_3=0.1,null=0)):
        '''
        example configuration is dictionary, see the `cfg` example above
        the proportions **do not have to sum to 1** because they will be normalized anyways
        '''
        default_keys = ['identity', 'equal_effects', 'singleton', 'simple_het_1', 'simple_het_2', 'simple_het_3', 'null']
        if set(list(cfg.keys())) != set(default_keys):
            raise ValueError(f"Input configuration has to contain all of these components: {default_keys}")
        sum_values = sum(cfg.values())
        for k in cfg:
            cfg[k] = cfg[k] / sum_values
        singleton = cfg.pop('singleton')
        for item in range(self.R):
            cfg[f'singleton_{item+1}'] = singleton / self.R
        self.pis.update(cfg)
        
    def apply_grid(self):
        def product(x,y):
            for item in y:
                yield x*item
        self.Us = dict(sum([[(f"{p}_{i+1}", g) for i, g in enumerate(product(self.Us[p], np.square(self.grid)))] for p in self.Us if p != 'null'], []) + \
                      [('null', self.Us['null'])])
        nG = len(self.grid)
        for k in list(self.pis.keys()):
            if k == 'null':
                continue
            for g in range(nG):
                if self.pis[k] > 0:
                    self.pis[f'{k}_{g+1}'] = self.pis[k] / nG
                else:
                    del self.Us[f'{k}_{g+1}']
            del self.pis[k]

    def get_effects(self):
        '''
        Generate B under multivariate normal mixture
        beta ~ \pi_0\delta_0 + \sum \pi_i N(0, U_i)
        '''
        self.coef = np.zeros((self.J, self.R))
        for j in range(self.J):
            # sample distribution
            dist_index = np.random.multinomial(1, list(self.pis.values()), size = 1).tolist()[0].index(1)
            name = list(self.pis.keys())[dist_index]
            self.coef[j,:] = np.random.multivariate_normal(self.mus['zeros'], self.Us[name], 1)
        
    def sparsify_effects(self, num_non_zero, top_only = False):
        '''
        top_only: only keep top `num_non_zero` effects
        '''
        if top_only:
            beta_max = np.amax(np.absolute(self.coef), axis = 1)
            big_beta_index = [i[0] for i in sorted(enumerate(beta_max), key = lambda x: x[1], reverse = True)]
            selected_index = big_beta_index[:min(len(big_beta_index), num_non_zero)]
        else:
            selected_index = np.random.choice(np.arange(self.J), size=num_non_zero)
        for j in range(self.J):
            if j not in selected_index:
                self.coef[j,:] = self.mus['zeros']
                
    def swap_top_effects(self, given_index):
        '''Set top effects to given indices
        One can specify index, or use the "top_index"
        generated by RegressionData.get_representative_features()
        '''
        given_index = np.array(given_index, dtype=int)
        nb = np.zeros(self.coef.shape)
        beta_max = np.amax(np.absolute(self.coef), axis = 1)
        big_beta_index = [i[0] for i in sorted(enumerate(beta_max), key = lambda x: x[1], reverse = True)]
        for idx in given_index:
            nb[idx,:] = self.coef[big_beta_index.pop(0),:]
        for idx in range(nb.shape[0]):
            if not idx in given_index:
                nb[idx,:] = self.coef[big_beta_index.pop(0),:]
        self.coef = nb
        
    def get_y(self, regression_data, pve, residual_corr, is_pve_per_variable = False):
        yhat = regression_data.X @ self.coef
        genetic_var = np.var(yhat,axis=0)
        if is_pve_per_variable:
            pve = [pve * np.count_nonzero(self.coef[:,i]) for i in range(self.R)]
            pve = np.array([x if x > 0 else 1 for x in pve])
        sigma = np.asarray(np.sqrt(genetic_var / pve - genetic_var)).flatten()
        # Some sigma might be zero if only one condition is regulated by a SNP.
        # We have to be able to set sigma to a non-zero (default) number for this case
        # a reasonable default might be the smallest of the rest of sigma?
        if 0 in sigma and np.any(sigma):
            sigma[sigma==0] = np.min(sigma[sigma!=0])
        if not np.any(sigma):
            sigma = sigma + 1
        sigma = np.diagflat(sigma)
        self.residual_variance = sigma @ residual_corr @ sigma
        return yhat + np.random.multivariate_normal(np.zeros(self.R), self.residual_variance, regression_data.X.shape[0])
    
    def get_prior(self):
        return dict([('xUlist', [self.Us[x] for x in self.Us if x != 'null']), 
                    ('pi', [self.pis[x] for x in self.Us if x != 'null']), 
                    ('null_weight', self.pis['null'])])

class MultivariateMixtureEqualEff(MultivariateMixture):
    '''In this simulation we assume all causal SNPs have fixed effects'''
    def __init__(self, dim):
        super().__init__(dim)

    def get_effects(self):
        '''
        Generate B under multivariate normal mixture
        beta ~ \pi_0\delta_0 + \sum \pi_i N(0, U_i)
        '''
        self.coef = np.zeros((self.J, self.R))
        dist_index = np.random.multinomial(1, list(self.pis.values()), size = 1).tolist()[0].index(1)
        name = list(self.pis.keys())[dist_index]
        effect = np.random.multivariate_normal(self.mus['zeros'], self.Us[name], 1)
        for j in range(self.J):
            # sample distribution
            self.coef[j,:] = effect
