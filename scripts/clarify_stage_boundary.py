"""Post hoc stage limits of the frozen score; see STAGE_BOUNDARY_CLARIFICATION_v1.md."""
import hashlib
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from run_analysis import ROOT, P, R, F
from extension_utils import extra_scores
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def stage_fit(label, c, s, levels, definitions):
    assert c.index.is_unique and s.index.is_unique
    ss = s.loc[c.index]
    assert ss[['Joint_program', 'Composite']].notna().all().all()
    X = pd.DataFrame({'const': np.ones(len(c))}, index=c.index)
    for lev in levels[1:]:
        X[lev] = c.boundary_stage.eq(lev).astype(float)
    X['Composite_z'] = (ss.Composite - ss.Composite.mean()) / ss.Composite.std(ddof=1)
    y = (ss.Joint_program - ss.Joint_program.mean()) / ss.Joint_program.std(ddof=1)
    assert np.linalg.matrix_rank(X) == X.shape[1]
    fit = sm.OLS(y, X).fit(cov_type='HC3', use_t=True)
    rows = []
    for name, vector in definitions:
        ct = fit.t_test(np.asarray(vector, dtype=float))
        ci = np.asarray(ct.conf_int()).ravel()
        # Algebra checked against explicit linear combinations of fitted stage means.
        effect = float(np.asarray(ct.effect).item())
        assert np.isclose(effect, np.dot(vector, fit.params))
        rows.append(dict(analysis=label, target='Joint_program', contrast=name,
                         n=len(c), residual_df=int(fit.df_resid), beta=effect,
                         ci_low=float(ci[0]), ci_high=float(ci[1]),
                         p=float(np.asarray(ct.pvalue).item()),
                         units='programme SD within this analysis sample set'))
    diag = dict(analysis=label, n=len(c), groups=c.boundary_stage.value_counts().to_dict(),
                design_rank=int(np.linalg.matrix_rank(X)), design_columns=X.columns.tolist(),
                residual_df=int(fit.df_resid), condition_number=float(np.linalg.cond(X)),
                is_independent_validation=False)
    points = c[['boundary_stage']].join(ss[['Joint_program', 'Composite']])
    points['analysis'] = label
    return rows, diag, points


def main():
    plan = ROOT / 'STAGE_BOUNDARY_CLARIFICATION_v1.md'
    assert plan.exists()
    plan_hash = hashlib.sha256(plan.read_bytes()).hexdigest()
    c = pd.read_csv(P / 'GSE63514_metadata.tsv', sep='\t', index_col=0)
    c = c.loc[c.whole_section.eq(0)].copy()
    c['boundary_stage'] = c.stage
    assert c.boundary_stage.value_counts().to_dict() == {'CIN3':40, 'Normal':24, 'CIN2':22, 'Cancer':20, 'CIN1':14}
    e = pd.read_csv(P / 'GSE63514_gene_log_expression.tsv.gz', sep='\t', index_col=0)
    s, _ = extra_scores(e)
    adjacent = [np.array(v, dtype=float) for v in (
        [0,1,0,0,0,0], [0,-1,1,0,0,0], [0,0,-1,1,0,0], [0,0,0,-1,1,0])]
    names = ['CIN1_vs_Normal','CIN2_vs_CIN1','CIN3_vs_CIN2','Cancer_vs_CIN3']
    defs = list(zip(names, adjacent)) + [('Cancer_CIN3_minus_'+name, adjacent[-1]-v)
                                        for name,v in zip(names[:3],adjacent[:3])]
    rows, diag, points = stage_fit('GSE63514_full_grades', c, s,
                                  ['Normal','CIN1','CIN2','CIN3','Cancer'], defs)
    diagnostics = [diag]
    all_points = [points]
    source_meta = pd.read_csv(P / 'GSE223804_metadata.tsv', sep='\t', index_col=0)
    matrix_to_gsm = source_meta['sample'].to_dict()
    c = source_meta.set_index('sample', drop=False)
    c['boundary_stage'] = np.select([c.stage.eq('HSIL'), c.title.str.contains('stage IA1|stage IA2', regex=True),
                                    c.stage.eq('Cancer')], ['HSIL','IA','IBplus'], default='Excluded')
    c = c.loc[c.boundary_stage.ne('Excluded')]
    assert c.boundary_stage.value_counts().to_dict() == {'IA':5, 'HSIL':4, 'IBplus':4}
    defs = [('IA_vs_HSIL',[0,1,0,0]),('IBplus_vs_HSIL',[0,0,1,0]),('IBplus_vs_IA',[0,-1,1,0])]
    e = pd.read_csv(P / 'GSE223804_log2_CPMplus1.tsv.gz', sep='\t', index_col=0)
    assert set(e.columns) == set(matrix_to_gsm)
    e = e.rename(columns=matrix_to_gsm)
    s, _ = extra_scores(e)
    for label, scores in [('GSE223804_source',s),
                          ('GSE223804_ARCHS4_same_patients',pd.read_csv(P/'GSE223804_ARCHS4_sum_scores.tsv', sep='\t', index_col=0))]:
        r, d, pts = stage_fit(label,c,scores,['HSIL','IA','IBplus'],defs)
        rows.extend(r); diagnostics.append(d); all_points.append(pts)
    result = pd.DataFrame(rows)
    assert len(result) == 13 and result.p.between(0,1).all()
    result['holm_13_stage_clarification'] = multipletests(result.p, method='holm')[1]
    result.to_csv(R/'stage_boundary_clarification.tsv', sep='\t', index=False)
    pd.concat(all_points).to_csv(R/'stage_boundary_sample_scores.tsv', sep='\t')
    (R/'stage_boundary_diagnostics.json').write_text(json.dumps({
        'plan_sha256':plan_hash,'diagnostics':diagnostics,'tests':13,
        'all_checks_passed':True,'interpretation':'Post hoc stage clarification; not new independent evidence. Different within-analysis SDs must not be read as directly comparable absolute effects.'},indent=2))
    fig, axs = plt.subplots(1,3,figsize=(13,3.8))
    for ax, pts, levels, title in zip(axs, all_points,
        [['Normal','CIN1','CIN2','CIN3','Cancer'],['HSIL','IA','IBplus'],['HSIL','IA','IBplus']],
        ['GSE63514 | laser capture','GSE223804 | source','GSE223804 | ARCHS4 (same patients)']):
        for i,lev in enumerate(levels):
            vals=pts.loc[pts.boundary_stage.eq(lev),'Joint_program'].to_numpy()
            ax.scatter(np.linspace(i-.13,i+.13,len(vals)),vals,s=22,alpha=.8,color=plt.cm.viridis((i+.5)/len(levels)))
            ax.plot([i-.18,i+.18],[np.median(vals)]*2,color='black',lw=1.4)
        ax.set(xticks=range(len(levels)),xticklabels=[f'{lev}\n(n={int(pts.boundary_stage.eq(lev).sum())})' for lev in levels],
               ylabel='Fixed signed mean fractional rank',title=title)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Cross-sectional stages: raw fixed scores (inference adjusts for proliferation)',fontsize=11)
    fig.tight_layout()
    fig.savefig(F/'20_stage_boundary_clarification.png',dpi=220)
    fig.savefig(F/'20_stage_boundary_clarification.pdf')
    plt.close(fig)
    print(result.to_string(index=False))
    print(json.dumps(diagnostics,indent=2))


if __name__ == '__main__':
    main()
