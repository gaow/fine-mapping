def summarize_LD(X, ld_input, ld_plot):
    data = RegressionData()
    data.X = X
    data.set_xcorr(ld_input)
    data.plot_xcorr(ld_plot)
    return data.get_representative_features()

def simulate_main(data, c, plot_prefix):
    '''
    data: $data
    top_idx: $top_eff
    n_signal: 3
    n_traits: 2
    eff_mode: mash_low_het
    swap_eff: True
    keep_ld: True
    tag: sim1
    @ALIAS: conf = Dict(!data, !eff_mode)
    $data: data
    '''
    reg = RegressionData()
    reg.X = data['X'].astype(float)
    if eff_mode == 'original':
        c['swap_eff'] = False
    if c['swap_eff'] and c['top_idx'] is None:
        raise ValueError(f'"top_idx" variable is not set by an upstream module')
    if eff_mode == 'mash_low_het':
        if c['n_traits'] < 2:
            raise ValueError(f'Cannot simulate {c["n_traits"]} under mode {eff_mode}')
        data['true_coef'], data['residual_variance'] = mash_low_het(data, reg, c)
    elif eff_mode == 'original':
        data['true_coef'] = original_y(data, reg, c)
        data['residual_variance'] = None
    elif eff_mode == 'simple_lm':
        data['true_coef'], data['residual_variance'] = simple_lm(data, reg, c)
    else:
        raise ValueError(f'Mode {eff_mode} is not implemented.')
    if c['center_data']:
        reg.center_data()
    data['X'] = reg.X
    data['Y'] = reg.Y
    data['allele_freq'] = (reg.x_mean / 2) if reg.x_mean is not None else (np.mean(reg.X, axis=0) / 2)
    data['allele_freq'] = data['allele_freq'].T
    if data['true_coef'] is not None:
        for j in range(data['true_coef'].shape[1]):
            plot_file = f'{plot_prefix}.{j+1}.png'
            reg.plot_property_vector(data['true_coef'][:,j], 
                                 [np.absolute(x)>0 for x in data['true_coef'][:,j]], 
                                 xz_cutoff = None, out = plot_file,
                                conf = {'title': f'Response {j+1}', 
                                        'ylabel': 'effect size', 'zlabel': ''})
    if data['Y'].shape[1] == 1:
        data['V'] = np.cov(data['Y'][:,0])
    else:
        data['V'] = np.cov(data['Y'], rowvar = False)
    return data
        
def original_y(data, reg, c):
    if 'y' in data:
        reg.Y = data.pop('y')
    elif isinstance(data['Y'], pd.DataFrame):
        reg.Y = np.vstack(data['Y'].values()).T
    else:
        reg.Y = data['Y']
    return None if 'true_coef' not in data else np.array(data['true_coef'])

def simple_lm(data, reg, c):
    if 'Z' in data:
        del data['Z']
    if 'y' in data:
        del data['y']
    eff = UnivariateMixture(reg.X.shape[1])
    eff.set_vanilla(c['amplitude'])
    Y = []
    coef = []
    sigma = []
    for i in range(c['n_traits']):
        eff.get_effects()
        eff.sparsify_effects(c['n_signal'])
        Y.append(eff.get_y(reg, pve = c['pve']))
        coef.append(eff.coef)
        sigma.append(eff.residual_variance)
    reg.Y = np.hstack(Y)
    return np.array(coef).T, np.array(sigma)