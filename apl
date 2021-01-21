#!/usr/bin/env python

import sys
import os
import argparse
import warnings

PATH = "git+https://github.com/lilyminium/mdanalysis@leaflet-methods#egg=mdanalysis&subdirectory=package"

try:
    import MDAnalysis as mda
except ImportError:
    raise ImportError(f'Install mdanalysis with `pip install -e "{PATH}"`')

try:
    from MDAnalysis.analysis.leaflets import LeafletFinder, AreaPerLipid
except ImportError:
    raise ImportError(f'Install mdanalysis with `pip install -e "{PATH}"`')


try:
    import pandas as pd
except ImportError:
    raise ImportError('Install pandas with `conda install pandas`')

try:
    import numpy as np
except ImportError:
    raise ImportError('Install numpy with `conda install numpy`')

# shut up
warnings.filterwarnings('ignore', message='Failed to guess the mass')

parser = argparse.ArgumentParser(description='Calculate area per lipid.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('topology', help='TPR or GRO or PDB file')
parser.add_argument('trajectory', help='XTC or TRR file')
parser.add_argument('--protein', default='protein', help='selection string for protein')
parser.add_argument('--lipids', default=['*'], 
                    nargs='+', 
                    help=('names of lipids to analyse, separated by space. '
                          'Accepts wildcards. e.g. --lipids P*PC CHOL'))
parser.add_argument('--headgroups', default=['PO4', 'GL1', 'ROH', 'GL2', 'AM1', 'AM2', 'NC3', 'GM*', 'C3'],
                    nargs='+', help=('names of headgroups, separated by space.'
                                     'Accepts wildcards. e.g. '
                                     '--headgroups ROH GL1 P*'))
parser.add_argument('--cutoff', default=50, type=float,
                    help='cutoff in ångström')
parser.add_argument('--leaflet_headgroups', default=['PO4', 'GL1', 'ROH', 'GL2', 'AM1', 'AM2', 'NC3', 'GM*', 'C3'],
                    nargs='+', help=('names of headgroups, separated by space.'
                                     'Accepts wildcards. e.g. '
                                     '--leaflet_headgroups ROH GL1 P*'))
parser.add_argument('--leaflet_method', default="spectralclustering", type=str,
                    help="leaflet grouping method")
parser.add_argument('--leaflet_cutoff', default=50, type=float,
                    help='leaflet cutoff in ångström')
parser.add_argument('--n_leaflets', default=2, type=int, help="number of leaflets")

# orientation
parser.add_argument('--angle_factor', default=1, type=float,
                    help="distance projection weighting (for orientation, spectralclustering method)")
parser.add_argument('--min_lipids', default=10, type=int,
                    help="minimum lipids per leaflet (for orientation method)")
parser.add_argument('--min_cosine', default=0.5, type=float,
                    help="cosine angle threshold between lipids(for orientation method)")
parser.add_argument('--max_neighbors', default=20, type=int,
                    help="maximum number of neighbors to group at once (for orientation method)")
parser.add_argument('--max_dist', default=20, type=float,
                    help="max radius to consider lipids in same leaflet (for orientation method)")
parser.add_argument('--relax_dist', default=10, type=float,
                    help="how much to relax max radius for lipid grouping (for orientation method)")

# spectral clustering
parser.add_argument('--delta', default=20, type=float,
                    help="distance weighting (for spectralclustering method)")
parser.add_argument('--angle_threshold', default=0.8, type=float,
                    help="how much to clip angles (for spectralclustering method)")

parser.add_argument('-b', '--begin', default=None, help='frame to begin at')
parser.add_argument('-e', '--end', default=None, help='frame to end at')
parser.add_argument('-skip', default=1, help='compute every nth frame')
parser.add_argument('--verbose', action='store_true', help='verbose')


if __name__ == '__main__':
    args = parser.parse_args(sys.argv[1:])
    u = mda.Universe(args.topology, args.trajectory)
    name, _ = os.path.splitext(os.path.basename(args.trajectory))

    
    headgroups = 'name ' + ' '.join(args.headgroups)
    SEL = 'name ' + ' '.join(args.leaflet_headgroups)

    lf = LeafletFinder(u, select=SEL, method=args.leaflet_method,
                       cutoff=args.leaflet_cutoff,
                       n_leaflets=args.n_leaflets,
                       angle_factor=args.angle_factor,
                       max_neighbors=args.max_neighbors,
                       min_cosine=args.min_cosine,
                       relax_dist=args.relax_dist,
                       min_lipids=args.min_lipids,
                       max_dist=args.max_dist, delta=args.delta,
                       angle_threshold=args.angle_threshold)

    apl = AreaPerLipid(u, leafletfinder=lf, select=headgroups,
                       verbose=args.verbose, select_other="not resname PW ION")
    print(f"Running analysis for {len(apl.residues)} residues. "
          "This could take a while. Perhaps put it in a "
          "tmux session and go grab a coffee for a week.")
    
    start = args.begin if not args.begin else int(args.begin)
    stop = args.end if not args.end else int(args.end)
    step = args.skip if not args.skip else int(args.skip)
    apl.run(start=start, stop=stop, step=step)
    name += f"_{apl.frames[0]:04d}-{apl.frames[-1]:04d}"

    leaflets = ["Extracellular", "Intracellular"]
    AREA = r"Area ($\AA^2$)"
    df_dct = {"Leaflet": [], "Lipid": [], AREA: []}
    for j, dct in enumerate(apl.areas_by_attr):
        for k, v in dct.items():
            n = len(v)
            df_dct["Leaflet"].extend([leaflets[j]]*n)
            df_dct["Lipid"].extend([k]*n)
            df_dct[AREA].extend(v)
    
    df = pd.DataFrame(df_dct)

    fn = f"{name}_apl_all.csv"
    df.to_csv(fn)
    print(f"Wrote {fn}")

    mean = df.groupby(["Leaflet", "Lipid"]).mean()
    mean.rename({AREA: r"Mean ($\AA^2$)"}, inplace=True)
    std = df.groupby(["Leaflet", "Lipid"]).std()
    mean[r"StDev ($\AA^2$)"] = std[AREA]

    print(mean)

    fn = f"{name}_apl_summary.csv"
    mean.to_csv(fn)
    print(f"Wrote {fn}")
    print('Done. ᕕ(⌐■_■)ᕗ ♪♬')

